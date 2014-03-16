import groovy.transform.CompileStatic;

import java.security.MessageDigest;


class Hash {

    public Hash() {
    }
    
    @CompileStatic
    static String sha1(String value) {
        def messageDigest = MessageDigest.getInstance("SHA1")
        messageDigest.update( value.getBytes() );
        return new BigInteger(1, messageDigest.digest()).toString(16).padLeft( 40, '0' )
    }
}
