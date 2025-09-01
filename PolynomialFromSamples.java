import org.json.JSONObject;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.math.BigInteger;
import java.util.*;

public class SecretFinder {

    // ---- minimal exact rational n/d (d>0, reduced) over BigInteger ----
    static final class Q {
        final BigInteger n, d;
        Q(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("denominator is zero");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            this.n = n.divide(g);
            this.d = d.divide(g);
        }
        static Q of(BigInteger n) { return new Q(n, BigInteger.ONE); }
        Q add(Q o) { return new Q(n.multiply(o.d).add(o.n.multiply(d)), d.multiply(o.d)); }
        Q mul(Q o) { return new Q(n.multiply(o.n), d.multiply(o.d)); }
        boolean isInt() { return d.equals(BigInteger.ONE); }
        @Override public String toString() { return isInt() ? n.toString() : (n + "/" + d); }
    }

    // ---- read whole file ----
    static JSONObject readJson(String path) throws Exception {
        String content = new String(Files.readAllBytes(Paths.get(path)));
        return new JSONObject(content);
    }

    // ---- decode y in base b (2..36) to BigInteger ----
    static BigInteger decode(String base, String value) {
        return new BigInteger(value, Integer.parseInt(base));
    }

    // ---- extract first k points (x=i, y decoded) in ascending i ----
    static void collectFirstK(JSONObject j, int k, List<Long> xs, List<BigInteger> ys) {
        int n = j.getJSONObject("keys").getInt("n");
        for (int i = 1; i <= n && xs.size() < k; i++) {
            String key = String.valueOf(i);
            if (!j.has(key)) continue; // skip missing indices just in case
            JSONObject pt = j.getJSONObject(key);
            xs.add((long)i);
            ys.add(decode(pt.getString("base"), pt.getString("value")));
        }
    }

    // ---- Lagrange interpolation at x = 0:  c = sum_i y_i * Π_{j≠i} (0 - x_j) / (x_i - x_j) ----
    static Q constantAtZero(List<Long> xs, List<BigInteger> ys) {
        int k = xs.size();
        Q c = Q.of(BigInteger.ZERO);
        for (int i = 0; i < k; i++) {
            BigInteger xi = BigInteger.valueOf(xs.get(i));
            BigInteger yi = ys.get(i);
            BigInteger num = BigInteger.ONE;   // Π (0 - xj)
            BigInteger den = BigInteger.ONE;   // Π (xi - xj)
            for (int j = 0; j < k; j++) {
                if (j == i) continue;
                BigInteger xj = BigInteger.valueOf(xs.get(j));
                num = num.multiply(xj.negate());
                den = den.multiply(xi.subtract(xj));
            }
            c = c.add(Q.of(yi).mul(new Q(num, den)));
        }
        return c;
    }

    // ---- compute secret c for a single JSON file ----
    static String solveFile(String path) throws Exception {
        JSONObject j = readJson(path);
        int k = j.getJSONObject("keys").getInt("k"); // degree m = k-1
        List<Long> xs = new ArrayList<>(k);
        List<BigInteger> ys = new ArrayList<>(k);
        collectFirstK(j, k, xs, ys);
        if (xs.size() < k) throw new IllegalArgumentException("Not enough points in " + path);
        Q c = constantAtZero(xs, ys);
        return c.toString(); // integer or "num/den"
    }

    public static void main(String[] args) {
        if (args.length == 0) {
            System.err.println("Usage: java SecretFinder <test1.json> [<test2.json> ...]");
            return;
        }
        for (String path : args) {
            try {
                String c = solveFile(path);
                System.out.println(c); // print just the secret, one per line (fits assignment output)
            } catch (Exception e) {
                System.err.println("Error in " + path + ": " + e.getMessage());
            }
        }
    }
}
