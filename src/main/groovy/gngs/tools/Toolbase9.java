package gngs.tools;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;

import sun.misc.Unsafe;

public class Toolbase9 {

    @SuppressWarnings({ "rawtypes", "unchecked" })
    public static void main(String[] args) {
        disableWarning();
        try {
            Class cls = Class.forName("gngs.tools." + args[0]);
            Method m = cls.getMethod("main", args.getClass());
            m.invoke(null, new Object[]{(Object[])Arrays.copyOfRange(args, 1, args.length)});
        } catch (ClassNotFoundException | NoSuchMethodException | SecurityException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            System.err.println("A problem was experienced executing the tool "+ args[0] + "\n\nThis may occur if you specified a tool that doesn't exist, or possibly if you are using an unusual JVM\n\n");
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
    
    /**
     * See https://stackoverflow.com/questions/46454995/how-to-hide-warning-illegal-reflective-access-in-java-9-without-jvm-argument
     */
    public static void disableWarning() {
        try {
            Field theUnsafe = Unsafe.class.getDeclaredField("theUnsafe");
            theUnsafe.setAccessible(true);
            Unsafe u = (Unsafe) theUnsafe.get(null);

            Class cls = Class.forName("jdk.internal.module.IllegalAccessLogger");
            Field logger = cls.getDeclaredField("logger");
            u.putObjectVolatile(cls, u.staticFieldOffset(logger), null);
        } catch (Exception e) {
            // ignore
        }
    }        

}
