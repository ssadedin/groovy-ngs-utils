package gngs

import groovy.json.JsonSlurper
import groovy.util.ProxyGenerator
import static org.junit.Assert.*

import org.junit.Test

import groovy.util.logging.Log

@Log
class WebserviceTest{

    // because HttpURLConnection is final in Java we cant metaClass it (or I couldnt, so here a solution for the problem)
    class MockHttpURLConnection extends HttpURLConnection {
        int callCount = 0
        int retries 
        ByteArrayOutputStream out = new ByteArrayOutputStream()

        MockHttpURLConnection(URL u, retries=0) { 
            super(u) 
            this.retries = retries
        }

        @Override void connect() {
            log.info "Mock connect()"
            callCount++
        }

        @Override int getResponseCode() {
            // response code can be queried multiple times, so we dont incremend the call count here
            if( callCount < retries){
                return 500
            }else{
                return 200
            }
        }

        @Override InputStream getInputStream() {
            return new ByteArrayInputStream('["Success after retry"]'.bytes)
        }

        @Override InputStream getErrorStream() {
            return new ByteArrayInputStream("Server error".bytes)
        }

        @Override String getContentType() {
            return "application/json"
        }

        @Override OutputStream getOutputStream() {
            return out
        }

        @Override void disconnect() {}
        @Override boolean usingProxy() { return false }
    }


    
    static {
        Utils.configureSimpleLogging()
        
    }
    
    @Test
    public void testRetryOnError() {
        def service = new WebService("http://mocked")
        def mock_conn = new MockHttpURLConnection(new URL("http://mocked"), 3)
        // Mock the configureConnection method to return the mock connection
        service.metaClass.configureConnection = { URL url, String method, Object data, Map headers ->
            return  mock_conn // Return the mocked connection
        }

        assert mock_conn.callCount == 0
        // add a new retry
        service.retry = new Retry()
        def before = System.currentTimeMillis()
        def r = service.get()
        def after = System.currentTimeMillis()

        // we should wait about 2 seconds
        assertEquals(2, (after-before)/1000, 0.1)
        // 3 tries with 1 success
        assert mock_conn.callCount == 3 
        assert r == ["Success after retry"]

    }

    @Test(expected = WebServiceException)
    public void testRetryButNotCode() {
        def service = new WebService("http://mocked")
        def mock_conn = new MockHttpURLConnection(new URL("http://mocked"), 3)
        // Mock the configureConnection method to return the mock connection
        service.metaClass.configureConnection = { URL url, String method, Object data, Map headers ->
            return  mock_conn // Return the mocked connection
        }

        service.retry = new Retry(retry_codes: [404])

        def r = service.get()
        
    }

    @Test(expected = WebServiceException)
    public void testNoRetry() {
        def service = new WebService("http://mocked")
        def mock_conn = new MockHttpURLConnection(new URL("http://mocked"), 3)
        // Mock the configureConnection method to return the mock connection
        service.metaClass.configureConnection = { URL url, String method, Object data, Map headers ->
            return  mock_conn // Return the mocked connection
        }

        def r = service.get()

    }

    @Test(expected = WebServiceException)
    public void testRunOutAttempts() {
        def service = new WebService("http://mocked")
        def mock_conn = new MockHttpURLConnection(new URL("http://mocked"), 5)
        // Mock the configureConnection method to return the mock connection
        service.metaClass.configureConnection = { URL url, String method, Object data, Map headers ->
            return  mock_conn // Return the mocked connection
        }

        service.retry = new Retry(total_retries: 1, back_off_factor: 0)
        def r = service.get()

    }


    @Test
    public void testRetryMaxSleep(){

        int max_sleep_seconds = 120
        long max_sleep_millis = max_sleep_seconds * 1000
        assertEquals max_sleep_millis, Retry.back_off_time(100, 1.0, max_sleep_seconds)

        
    }

    @Test
    public void testOAuth2TokenAcquisition() {
        def service = new WebService("http://api.example.com")
        
        // Mock the token endpoint
        def tokenCallCount = 0
        
        // Mock the API endpoint
        def apiCallCount = 0
        def mockApiConnection = new MockHttpURLConnection(new URL("http://api.example.com")) {
            @Override
            InputStream getInputStream() {
                apiCallCount++
                // Verify the bearer token was set
                assert getRequestProperty('Authorization') == 'Bearer test-token-123'
                return new ByteArrayInputStream('["API Success"]'.bytes)
            }
        }
        
        // Create OAuth2 credentials with mocked token fetcher
        def oauth2Creds = new OAuth2ClientCredentials(
            clientId: 'test-client',
            clientSecret: 'test-secret',
            tokenUrl: 'http://auth.example.com/token'
        )
        
        // Mock the token fetcher to avoid network calls
        oauth2Creds.tokenFetcher = {
            tokenCallCount++
            return [access_token: 'test-token-123', expires_in: 3600]
        }
        
        service.oauth2Credentials = oauth2Creds
        
        // Mock configureConnection to use our mock API connection
        service.metaClass.configureConnection = { URL url, String method, Object data, Map headers ->
            oauth2Creds.configure(mockApiConnection, url, method, data, headers)
            return mockApiConnection
        }
        
        def result = service.get([:])
        
        // Verify token was acquired once
        assert tokenCallCount == 1
        // Verify API was called once
        assert apiCallCount == 1
        assert result == ["API Success"]
    }

    @Test
    public void testOAuth2TokenRefresh() {
        def service = new WebService("http://api.example.com")
        
        def tokenCallCount = 0
        
        def apiCallCount = 0
        def capturedTokens = []
        def mockApiConnection = new MockHttpURLConnection(new URL("http://api.example.com")) {
            @Override
            InputStream getInputStream() {
                apiCallCount++
                def authHeader = getRequestProperty('Authorization')
                capturedTokens << authHeader
                return new ByteArrayInputStream('["API Success"]'.bytes)
            }
        }
        
        def oauth2Creds = new OAuth2ClientCredentials(
            clientId: 'test-client',
            clientSecret: 'test-secret',
            tokenUrl: 'http://auth.example.com/token'
        )
        
        oauth2Creds.tokenFetcher = {
            tokenCallCount++
            return [access_token: "test-token-${tokenCallCount}", expires_in: 1]
        }
        
        service.oauth2Credentials = oauth2Creds
        
        service.metaClass.configureConnection = { URL url, String method, Object data, Map headers ->
            oauth2Creds.configure(mockApiConnection, url, method, data, headers)
            return mockApiConnection
        }
        
        // First call - should acquire token
        service.get([:])
        assert tokenCallCount == 1
        assert capturedTokens[0] == 'Bearer test-token-1'
        
        // Second call - token should be expired (expires_in=1 means already expired with 60s buffer)
        // Force expiry by setting expiresAt to past
        oauth2Creds.@expiresAt = System.currentTimeMillis() - 1000
        
        service.get([:])
        assert tokenCallCount == 2
        assert capturedTokens[1] == 'Bearer test-token-2'
    }

    @Test(expected = WebServiceException)
    public void testOAuth2TokenAcquisitionFailure() {
        def service = new WebService("http://api.example.com")
        
        def oauth2Creds = new OAuth2ClientCredentials(
            clientId: 'bad-client',
            clientSecret: 'bad-secret',
            tokenUrl: 'http://auth.example.com/token'
        )
        
        oauth2Creds.tokenFetcher = {
            throw new WebServiceException(
                "Failed to obtain OAuth2 token from ${oauth2Creds.tokenUrl}: {\"error\":\"invalid_client\"}",
                401,
                'Unauthorized',
                '{"error":"invalid_client"}'
            )
        }
        
        service.oauth2Credentials = oauth2Creds
        
        // This should trigger token acquisition which will fail
        service.metaClass.configureConnection = { URL url, String method, Object data, Map headers ->
            oauth2Creds.configure(new MockHttpURLConnection(url), url, method, data, headers)
            return new MockHttpURLConnection(url)
        }
        
        // Should throw WebServiceException
        service.get([:])
    }

    @Test
    public void testOAuth2WithScope() {
        def oauth2Creds = new OAuth2ClientCredentials(
            clientId: 'test-client',
            clientSecret: 'test-secret',
            tokenUrl: 'http://auth.example.com/token',
            scope: 'read write'
        )
        
        def capturedScope = null
        oauth2Creds.tokenFetcher = {
            // Capture the scope to verify it would be used
            capturedScope = oauth2Creds.scope
            return [access_token: 'scoped-token', expires_in: 3600]
        }
        
        // Trigger token acquisition
        oauth2Creds.ensureValidToken()
        
        // Verify scope was captured
        assert capturedScope == 'read write'
        assert oauth2Creds.@accessToken == 'scoped-token'
    }
}
