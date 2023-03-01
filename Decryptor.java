package ckks;

public class Decryptor {
    private SecretKey secretKey;

    public Decryptor(SecretKey secretKey) {
        this.secretKey = secretKey;
    }

    public void decrypt(Ciphertext src, Plaintext res) {
        if (src.isEmpty())
            throw new IllegalArgumentException("Ciphertext is empty.");

        Polynomial s = secretKey.getS();

        Polynomial m;

        if (src.isLinear()) {
            Polynomial b = src.getB();
            Polynomial a = src.getA();

            m = b.add(a.mult(s));
        } else {
            Polynomial d0 = src.getD0();
            Polynomial d1 = src.getD1();
            Polynomial d2 = src.getD2();

            m = d0.add(d1.mult(s)).add(d2.mult(s.square()));
        }

        res.init(m, src.getScale(), src.getLevel());
    }
}
