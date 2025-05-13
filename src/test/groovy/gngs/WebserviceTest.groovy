package gngs

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
}