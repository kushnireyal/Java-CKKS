package ckks;

public class SecretKey {
    private Polynomial s;

    public SecretKey() {
    }

    public SecretKey(Polynomial s) {
        this.s = s;
    }

    public Polynomial getS() {
        return s;
    }

    public long[] serialize() {
        return s.serialize();
    }

    public void deserialize(Context context, long[] serialization) {
        s = new Polynomial();
        s.deserialize(context, serialization);
    }
}
