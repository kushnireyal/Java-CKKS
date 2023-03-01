package ckks;

public class Ciphertext {
    private Context context;

    private Polynomial b, a;

    private Polynomial d0, d1, d2;

    private double scale;

    private int level;

    public Ciphertext(Context context) {
        this.context = context;
    }

    public Ciphertext(Ciphertext src) {
        this.context = src.context;

        if (src.isEmpty())
            return;

        if (src.isLinear()) {
            this.b = new Polynomial(src.b);
            this.a = new Polynomial(src.a);
        } else {
            this.d0 = new Polynomial(src.d0);
            this.d1 = new Polynomial(src.d1);
            this.d2 = new Polynomial(src.d2);
        }

        this.scale = src.scale;

        this.level = src.level;
    }

    public void init(Polynomial b, Polynomial a, int level,
            double scale) {
        this.b = b;
        this.a = a;
        this.level = level;
        this.scale = scale;
    }

    public void afterMult(Polynomial d0, Polynomial d1, Polynomial d2,
            double scale) {
        this.b = null;
        this.a = null;
        this.d0 = d0;
        this.d1 = d1;
        this.d2 = d2;
        this.scale = scale;
    }

    public void afterRelinearization(Polynomial b, Polynomial a) {
        this.d0 = null;
        this.d1 = null;
        this.d2 = null;
        this.b = b;
        this.a = a;
    }

    public void afterRescale(Polynomial b, Polynomial a, double scale,
            int level) {
        this.d0 = null;
        this.d1 = null;
        this.d2 = null;
        this.b = b;
        this.a = a;
        this.scale = scale;
        this.level = level;
    }

    public Polynomial getB() {
        return b;
    }

    public Polynomial getA() {
        return a;
    }

    public Polynomial getD0() {
        return d0;
    }

    public Polynomial getD1() {
        return d1;
    }

    public Polynomial getD2() {
        return d2;
    }

    public boolean isEmpty() {
        return b == null && d0 == null;
    }

    public boolean isLinear() {
        return b != null;
    }

    public double getScale() {
        return scale;
    }

    public void setScale(double scale) {
        this.scale = scale;
    }

    public int getLevel() {
        return level;
    }

    public void setLevel(int level) {
        this.level = level;
    }

    public void debugPrint() {
        if (isLinear()) {
            System.out.println("Ciphertext is linear");
            System.out.println("b:");
            b.debugPrint();
            System.out.println("a:");
            a.debugPrint();
        } else {
            System.out.println("Ciphertext is not linear");
            System.out.println("d0:");
            d0.debugPrint();
            System.out.println("d1:");
            d1.debugPrint();
            System.out.println("d2:");
            d2.debugPrint();
        }
    }

    public long[] serialize() {
        if (d0 != null)
            throw new IllegalArgumentException("Ciphertext::serialize supported only for linear ciphertext.");

        long[][] bCrt = b.getCrt();
        long[][] aCrt = a.getCrt();

        long[] res = new long[(2 + bCrt.length * bCrt[0].length) + (2 + aCrt.length * aCrt[0].length) + 2];
        int idx = 0;

        res[idx++] = bCrt.length;
        res[idx++] = bCrt[0].length;
        for (int i = 0; i < bCrt.length; i++)
            for (int j = 0; j < bCrt[0].length; j++)
                res[idx++] = bCrt[i][j];

        res[idx++] = aCrt.length;
        res[idx++] = aCrt[0].length;
        for (int i = 0; i < aCrt.length; i++)
            for (int j = 0; j < aCrt[0].length; j++)
                res[idx++] = aCrt[i][j];

        res[idx++] = (long) scale;
        res[idx] = level;

        return res;
    }

    public void deserialize(long[] serialization) {
        int idx = 0;

        long[][] bCrt = new long[(int) serialization[idx++]][(int) serialization[idx++]];
        for (int i = 0; i < bCrt.length; i++)
            for (int j = 0; j < bCrt[0].length; j++)
                bCrt[i][j] = serialization[idx++];

        long[][] aCrt = new long[(int) serialization[idx++]][(int) serialization[idx++]];
        for (int i = 0; i < aCrt.length; i++)
            for (int j = 0; j < aCrt[0].length; j++)
                aCrt[i][j] = serialization[idx++];

        b = new Polynomial(context, bCrt);
        a = new Polynomial(context, aCrt);

        d0 = null;
        d1 = null;
        d2 = null;

        scale = serialization[idx++];
        level = (int) serialization[idx];
    }
}
