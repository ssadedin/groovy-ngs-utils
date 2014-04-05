
class Time {
    
    static Hashtable<String,Stats> timings = new Hashtable()

    public Time() {
    }
    
    static void time(String key, Closure c) {
        long startTimeMs = System.currentTimeMillis()
        c()
        long endTimeMs = System.currentTimeMillis()
        
        Stats s = timings[key]
        if(s == null) {
            s = new Stats()
            timings[key] = s
        }
        
        timings[key].addValue(endTimeMs - startTimeMs)
    }
    
    static void dump() {
        
        // longest key to determine column width
        int width = timings.keySet().max { it.size() }.size() + 2
        
        int dataWidth = 12
        println "=" * 100
        println "Key".padRight(width) + "Calls".padRight(dataWidth) + "Mean Time".padRight(dataWidth)
        timings.keySet().each { key ->
            Stats s = timings[key]
            println key.padRight(width) + s.getN().toString().padRight(dataWidth) + s.mean.toString().padRight(dataWidth)
        }
        
    }

}
