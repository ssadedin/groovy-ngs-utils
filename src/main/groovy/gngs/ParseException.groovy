package gngs

class ParseException extends Exception {
    
    Integer line

    public ParseException() { super(); }

    public ParseException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) { super(message, cause, enableSuppression, writableStackTrace); }

    public ParseException(String message, Throwable cause) {
        super(message, cause);
    }

    public ParseException(String message, Integer line) { super("Failure at line " + line + ": " + message); this.line = line; }

    public ParseException(String message, Integer line, Throwable cause) { super("Failure at line " + line + ": " + message, cause); this.line = line; }

    public ParseException(Throwable cause, Integer line) { super("Failure at line $line: " + cause.getMessage(), cause); this.line = line; }

    public ParseException(Throwable cause) { super(cause); }
    
}
