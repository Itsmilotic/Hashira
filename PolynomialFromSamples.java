import java.io.*;
import java.math.BigInteger;
import java.util.*;
import org.json.JSONObject;

public class PolynomialFromSamples {
    // ---- Minimal BigInteger Rational ----
    static final class Q {
        BigInteger n, d; // always d>0, gcd=1git
        Q(BigInteger n, BigInteger d){ if(d.signum()==0) throw new ArithmeticException("den=0");
            if(d.signum()<0){ n=n.negate(); d=d.negate(); }
            BigInteger g = n.gcd(d); this.n = n.divide(g); this.d = d.divide(g);
        }
        Q(long n){ this(BigInteger.valueOf(n), BigInteger.ONE); }
        static Q of(BigInteger n){ return new Q(n, BigInteger.ONE); }
        Q add(Q o){ return new Q(n.multiply(o.d).add(o.n.multiply(d)), d.multiply(o.d)); }
        Q sub(Q o){ return new Q(n.multiply(o.d).subtract(o.n.multiply(d)), d.multiply(o.d)); }
        Q mul(Q o){ return new Q(n.multiply(o.n), d.multiply(o.d)); }
        Q div(Q o){ return new Q(n.multiply(o.d), d.multiply(o.n)); }
        boolean isZero(){ return n.signum()==0; }
        public String toString(){ return d.equals(BigInteger.ONE) ? n.toString() : (n + "/" + d); }
    }

    // Multiply polynomial A by (x - a): coeffs are ascending degree
    static List<Q> mulByXMinusA(List<Q> A, BigInteger a){
        int m = A.size();
        List<Q> C = new ArrayList<>(Collections.nCopies(m+1, new Q(0)));
        for(int i=0;i<m;++i){
            // C[i] += A[i]*(-a)
            C.set(i, C.get(i).add(A.get(i).mul(new Q(a.negate(), BigInteger.ONE))));
            // C[i+1] += A[i]
            C.set(i+1, C.get(i+1).add(A.get(i)));
        }
        return C;
    }
    static List<Q> addPoly(List<Q> A, List<Q> B){
        int m=Math.max(A.size(),B.size());
        List<Q> C=new ArrayList<>(Collections.nCopies(m,new Q(0)));
        for(int i=0;i<m;++i){
            Q ai = i<A.size()?A.get(i):new Q(0);
            Q bi = i<B.size()?B.get(i):new Q(0);
            C.set(i, ai.add(bi));
        }
        return C;
    }
    static List<Q> scalePoly(List<Q> A, Q s){
        List<Q> B=new ArrayList<>(A.size());
        for(Q c:A) B.add(c.mul(s));
        return B;
    }

    static BigInteger parseInBase(String s, int base){
        return new BigInteger(s, base); // handles 2..36, case-insensitive
    }

    public static void main(String[] args) throws Exception {
        // Read JSON from STDIN (fallback to sample if none given)
        String input;
        try(Scanner sc=new Scanner(System.in).useDelimiter("\\A")){
            input = sc.hasNext()? sc.next() : "";
        }
        if(input.isEmpty()){
            input = "{ \"keys\":{\"n\":4,\"k\":3},"
                  + "\"1\":{\"base\":\"10\",\"value\":\"4\"},"
                  + "\"2\":{\"base\":\"2\",\"value\":\"111\"},"
                  + "\"3\":{\"base\":\"10\",\"value\":\"12\"},"
                  + "\"6\":{\"base\":\"4\",\"value\":\"213\"}}";
        }
        JSONObject obj = new JSONObject(input);
        int n = obj.getJSONObject("keys").getInt("n");
        int k = obj.getJSONObject("keys").getInt("k"); // degree m=k-1

        // Build points (x=i, y parsed). Use i=1..n to keep x in order.
        List<Integer> xs = new ArrayList<>();
        List<BigInteger> ys = new ArrayList<>();
        for(int i=1;i<=n;i++){
            String key = Integer.toString(i);
            if(!obj.has(key)) continue;
            JSONObject pt = obj.getJSONObject(key);
            int base = Integer.parseInt(pt.getString("base"));
            String val = pt.getString("value");
            xs.add(i);
            ys.add(parseInBase(val, base));
        }
        if(xs.size() < k){ System.err.println("Not enough points."); return; }

        // Use first k points (interpolate degree k-1 polynomial)
        xs = xs.subList(0,k);
        ys = ys.subList(0,k);

        // Lagrange interpolation (exact rationals)
        List<Q> res = new ArrayList<>(); res.add(new Q(0));
        for(int i=0;i<k;i++){
            // basis = Π_{j≠i} (x - xj)
            List<Q> basis = new ArrayList<>(); basis.add(new Q(1));
            BigInteger denom = BigInteger.ONE;
            int xi = xs.get(i);
            for(int j=0;j<k;j++){
                if(j==i) continue;
                int xj = xs.get(j);
                basis = mulByXMinusA(basis, BigInteger.valueOf(xj));
                denom = denom.multiply(BigInteger.valueOf(xi - xj));
            }
            Q scale = new Q(ys.get(i), denom); // yi / Π(xi - xj)
            res = addPoly(res, scalePoly(basis, scale));
        }
        // trim trailing zeros
        while(res.size()>1 && res.get(res.size()-1).isZero()) res.remove(res.size()-1);

        // Output
        int degree = res.size()-1;
        System.out.println("{");
        System.out.println("  \"degree\": " + degree + ",");
        System.out.print("  \"coefficients_high_to_low\": [");
        for(int i=degree;i>=0;i--){
            System.out.print("\""+res.get(i).toString()+"\"");
            if(i>0) System.out.print(", ");
        }
        System.out.println("]\n}");
    }
}
