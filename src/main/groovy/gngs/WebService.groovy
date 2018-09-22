/*
 *  Groovy NGS Utils - Utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

package gngs

import org.yaml.snakeyaml.Yaml

import com.github.scribejava.core.builder.ServiceBuilder
import com.github.scribejava.core.builder.api.DefaultApi10a
import com.github.scribejava.core.model.OAuth1AccessToken
import com.github.scribejava.core.model.OAuthRequest
import com.github.scribejava.core.model.Response
import com.github.scribejava.core.model.Verb
import com.github.scribejava.core.oauth.OAuth10aService

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.util.logging.Log
import java.text.ParseException

/**
 * With 2 legged OAuth there is no need for acquiring an access token or request token; these
 * are presumed known. Therefore we can implement these methods with dummy bodyies.
 * 
 * @author Simon Sadedin
 */
class TwoLeggedOAuthAPI extends DefaultApi10a {

    @Override
    public String getAccessTokenEndpoint() {
        return null;
    }

    @Override
    protected String getAuthorizationBaseUrl() {
        return null;
    }

    @Override
    public String getRequestTokenEndpoint() {
        return null;
    }
}

class WSCredentials  {
    String apiKey
    String apiSecret
}

class WebServiceException extends Exception {
    public WebServiceException(String message, int code, String reason, String body) {
        super(message);
        this.code = code;
        this.reason = reason;
        this.body = body;
    }
    
    public WebServiceException(Throwable t) {
        super('Unknown error occurred executing request', t)
    }
    
    int code
    String reason
    String body
}

/**
 * Base class for APIs accessing Web services based on JSON data
 *
 * <p>
 * Usage:
 * <pre>
 * WebService w = new WebService("http://localhost:8000","some/url/path")
 * println "Result of call is: " + w.get(foo:10,bar:"tree")
 * </pre>
 *
 * @author simon.sadedin
 */
@Log
class WebService {
    
    String endPoint
    
    String api
    
    OAuth1AccessToken credentials
    
    WSCredentials apiCredentials
    
    String secret
    
    boolean verbose = true
    
    boolean dump = false
    
    boolean autoSlash = false
    
    int maxBodyDumpSize = 1024
    
    String credentialsPath = null
    
    /**
     * Create a web service to access the given api at the given end point
     *
     * @param endPoint  root URL of the webservice implementation (eg: http://localhost:8000)
     * @param api       api path, eg: 'xip/sample/16W0000001'
     */
    public WebService(String endPoint, String api=null) {
        super();
        this.endPoint = endPoint;
        if(api != null)
            this.api = api.startsWith('/') ? api.substring(1) : api;
        else
            this.api = ''
    }
    
    Object post(Object data) {
        request([:], 'POST', data)
    }
  
    Object post(Map params, Object data) {
        request(params, 'POST', data)
    }
    
    Object get(Map params) {
        request(params, 'GET', null)
    }
    
    WebService path(String subPath) {
        WebService child = new WebService(endPoint, api + '/' + subPath)
        child.verbose = this.verbose
        child.dump = this.dump
        child.credentials = this.credentials
        child.apiCredentials = this.apiCredentials
        child.autoSlash = this.autoSlash
        child.credentialsPath = this.credentialsPath
        return child
    }
    
    WebService div(String subPath) {
        path(subPath)
    }
    
    /**
     * Execute the given request using the given method, sending the given data as POST
     * data, if provided.
     *
     * @return an object representing the result.
     */
    Object request(Map params, String method, Object data) {
        
       String payload = data != null ? JsonOutput.prettyPrint(JsonOutput.toJson(data)) : null
       
       URL url = encodeURL(params, payload)
       
       try {
           if(credentials) 
               return this.executeOAuthRequest(params, url, method, payload)
           else {
               HttpURLConnection connection = configureConnection(url, method, data)
               return executeRequest(connection, payload)
           }
       }
       catch(WebServiceException e) {
           log.severe "Request to URL $url failed"
           if(dump && e.body) {
               File responseFile = new File(api.replaceAll('/','_'))
               log.severe("Dumping response body to $responseFile")
               responseFile.text = e.body
           }
           
           if(e.body) {
               if(e.body.size()<maxBodyDumpSize)
                   log.severe("Response content:\n" + e.body)
               else
                   log.severe("Body contains ${e.body.size()} chars, too large to dump")
           }
               
           throw e
       }
       catch(Exception e) {
           log.severe "Request to URL $url failed. Payload: \n$payload"
           throw new WebServiceException(e)
       }
    }
    
    Object executeOAuthRequest(Map params, URL url, String method, String data) {
        
        Verb verb = Verb.GET;
        
        switch(method) {
            case "GET": verb = Verb.GET; break
            case "POST": verb = Verb.POST; break
            case "PUT": verb = Verb.PUT; break
            case "DELETE": verb = Verb.DELETE; break
        }
        
        OAuth10aService service = buildOAuthService()

        log.info "Sending URL $url"
        OAuthRequest request = new OAuthRequest(verb, url.toString());
        request.addHeader('Accept', 'application/json')
        service.signRequest(credentials, request)
        if(data != null) {
            request.addHeader('Content-Type', 'application/json')
            request.setPayload(data)
        }
        Response response = service.execute(request)
        if(response.getCode() >= 400) {
            WebServiceException e = 
                new WebServiceException("Request to $url failed with status code $response.code (response=${response.body?.take(80)}...)", response.code, response.getMessage(), response.body)
            throw e
        }
        String responseText = response.body
        return convertResponse(response.getHeader('Content-Type'), responseText)
    }
    
    OAuth10aService buildOAuthService() {
        
        loadCredentials()

        new ServiceBuilder(this.apiCredentials.apiKey)
            .apiSecret(this.apiCredentials.apiSecret)
            .build(new TwoLeggedOAuthAPI()) 
    }

    public void loadCredentials() {
        if(this.apiCredentials != null) 
            return 
            
        File credsFile 
        List<File> credsFiles = []
        if(credentialsPath)  {
            credsFiles =[new File(credentialsPath), new File(System.properties['user.home'],credentialsPath) ]
            credsFile = credsFiles.find { it.exists() }
        }
        if(!credsFile)
            throw new IllegalStateException("API Credentials were not set explicitly and no credentials file could be found in ${credsFiles.join(', ')}")

        def yaml = new Yaml().load(credsFile.text)
        
        if(!(yaml instanceof Map)) 
            throw new ParseException("Bad format of credentials file. Please use YAML style syntax to define key value pairs")
        
        log.info "Loaded credentials from credentials file"
        
        if(yaml.apiKey) {
            this.apiCredentials = new WSCredentials(apiKey:yaml.apiKey, apiSecret: yaml.apiSecret)
        }
        
        if(yaml.accessToken) {
            this.credentials = new OAuth1AccessToken(yaml.accessToken, yaml.accessSecret)
        }
    }
    
    /**
     * Set up and return a connection configured to execute the given request
     */
    HttpURLConnection configureConnection(URL url, String method, Object data) {
        HttpURLConnection connection = url.openConnection()
        connection.with {
            doOutput = (data != null)
            useCaches = false
            setRequestProperty('Content-Type','application/json')
            requestMethod = method
        }
        return connection
    }
    
    /**
     * @return  an object whose type depends on the content type returned by the remote call.
     *          If the content type is JSON, the JSON is decoded and the object corresponds to that of the
     *          decoded content. Otherwise, the object is a string.
     */
    Object executeRequest(HttpURLConnection connection, String body) {
        connection.connect()
        if(body != null) {
            connection.outputStream.withWriter { writer ->
              writer << body
            }
        }
            
        if(verbose)
            log.info "Sent $connection.requestMethod to URL $connection.url"
                    
        int code = connection.getResponseCode()
        if(verbose)
            log.info("Received response code $code from server")
                
        String responseText = (code < 400) ? connection.inputStream.text : connection.errorStream.text
        if(code >= 400) {
            WebServiceException e = new WebServiceException("Request to $connection.url failed with status code $code (response=${responseText?.take(80)}...)", code, connection.responseMessage, responseText)
//            if(body != null)
//                e.body = body
            throw e
        }
                
        return convertResponse(connection.getContentType(), responseText)
    }
    
    /**
     * Convert the response to an appropriate type based on the response content-type
     * header.
     *
     * @param connection    Connection used for the transaction
     * @param responseText  response output
     *
     * @return  an object that is corresponds to the parsed form of the response, with format
     *          indicated by the content-type header.
     */
    Object convertResponse(String contentType, String responseText) {
        if(contentType == 'application/json') {
            try {
                if(responseText.isEmpty())
                    return null
                return new JsonSlurper().parseText(responseText)
            }
            catch(Exception e) {
                throw new WebServiceException("Failed to parse response as JSON: ${responseText?.take(80)}...", connection.responseCode, connection.responseMessage, responseText)
            }
        }
        else {
            return responseText
        }
    }
    
    URL encodeURL(Map params, String payload) {
        String url = "$endPoint/$api"
        if(credentials) {
            if(payload) {
                if(params == null)
                    params = [:]
                params.sig = Hash.sha1(payload, 'UTF-8')
            } 
        }
        
        if(this.autoSlash) {
            String lastElement = url.tokenize('/')[-1]
            if(!lastElement.contains('.') && !url.endsWith('/')) {
                url = url + '/'
            }
        }
        
        if(params) {
             url += "?" + params.collect { key, value ->
                               URLEncoder.encode(key) + '=' + URLEncoder.encode(String.valueOf(value))
                          }.join('&')
        }
        new URL(url)
    }
}

