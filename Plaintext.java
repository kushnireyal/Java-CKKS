package ckks;

public class Plaintext {
    private Context context;

    private Polynomial m;

    private double scale;

    private int level;

    public Plaintext(Context context) {
        this.context = context;
    }

    public Plaintext(Plaintext src) {
        this.context = src.context;

        if (src.isEmpty())
            return;

        this.m = new Polynomial(src.m);

        this.scale = src.scale;

        this.level = src.level;
    }

    public void init(Polynomial m, int level) {
        init(m, context.defaultScale, level);
    }

    public void init(Polynomial m, double scale, int level) {
        this.m = m;
        this.level = level;
        this.scale = scale;
    }

    public Polynomial getM() {
        return m;
    }

    public boolean isEmpty() {
        return m == null;
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
        System.out.println("m:");
        m.debugPrint();
    }
}
