package ckks;

public class PublicKeys {
    private Polynomial b, a;

    long[][] relinKeyB, relinKeyA;

    public PublicKeys() {
    }

    public PublicKeys(Polynomial b, Polynomial a) {
        this.b = b;
        this.a = a;
    }

    public void setRelinKeys(long[][] relinKeyB,
            long[][] relinKeyA) {
        this.relinKeyB = relinKeyB;
        this.relinKeyA = relinKeyA;
    }

    public Polynomial getB() {
        return b;
    }

    public Polynomial getA() {
        return a;
    }

    public long[][] getRelinKeyB() {
        return relinKeyB;
    }

    public long[][] getRelinKeyA() {
        return relinKeyA;
    }

    public long[] serialize() {
        long[] bSerialized = b.serialize();
        long[] aSerialized = a.serialize();

        long[] res = new long[1 + bSerialized.length + 1 + aSerialized.length + 2
                + relinKeyA.length * relinKeyA[0].length + 2 + relinKeyB.length * relinKeyB[0].length];

        int idx = 0;

        res[idx++] = bSerialized.length;
        for (int i = 0; i < bSerialized.length; i++)
            res[idx++] = bSerialized[i];

        res[idx++] = aSerialized.length;
        for (int i = 0; i < aSerialized.length; i++)
            res[idx++] = aSerialized[i];

        res[idx++] = relinKeyB.length;
        res[idx++] = relinKeyB[0].length;
        for (int i = 0; i < relinKeyB.length; i++)
            for (int j = 0; j < relinKeyB[0].length; j++)
                res[idx++] = relinKeyB[i][j];

        res[idx++] = relinKeyA.length;
        res[idx++] = relinKeyA[0].length;
        for (int i = 0; i < relinKeyA.length; i++)
            for (int j = 0; j < relinKeyA[0].length; j++)
                res[idx++] = relinKeyA[i][j];

        return res;
    }

    public void deserialize(Context context, long[] serialization) {
        int idx = 0;

        long[] bSerialized = new long[(int) serialization[idx++]];
        for (int i = 0; i < bSerialized.length; i++)
            bSerialized[i] = serialization[idx++];
        b = new Polynomial();
        b.deserialize(context, bSerialized);

        long[] aSerialized = new long[(int) serialization[idx++]];
        for (int i = 0; i < aSerialized.length; i++)
            aSerialized[i] = serialization[idx++];
        a = new Polynomial();
        a.deserialize(context, aSerialized);

        relinKeyB = new long[(int) serialization[idx++]][(int) serialization[idx++]];
        for (int i = 0; i < relinKeyB.length; i++)
            for (int j = 0; j < relinKeyB[0].length; j++)
                relinKeyB[i][j] = serialization[idx++];

        relinKeyA = new long[(int) serialization[idx++]][(int) serialization[idx++]];
        for (int i = 0; i < relinKeyA.length; i++)
            for (int j = 0; j < relinKeyA[0].length; j++)
                relinKeyA[i][j] = serialization[idx++];
    }
}
