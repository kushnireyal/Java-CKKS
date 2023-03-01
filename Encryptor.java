package ckks;

public class Encryptor {
    final boolean debug = false;

    private Context context;
    PublicKeys publicKey;

    public Encryptor(Context context, PublicKeys publicKey) {
        this.context = context;
        this.publicKey = publicKey;
    }

    public void encrypt(Plaintext src, Ciphertext res) {
        Polynomial v = KeyGenerator.ternaryDist(context, 0.5);
        Polynomial e1 = KeyGenerator.ternaryDist(context, 0.5);
        Polynomial e2 = KeyGenerator.ternaryDist(context, 0.5);

        if (debug) {
            System.out.println("Encryptor.encrypt");
            System.out.println("v");
            v.debugPrint();
            System.out.println("");
            System.out.println("e1");
            e1.debugPrint();
            System.out.println("");
            System.out.println("e2");
            e2.debugPrint();
        }

        Polynomial m = src.getM();
        Polynomial bKey = publicKey.getB();
        Polynomial aKey = publicKey.getA();

        Polynomial b = v.mult(bKey).add(m).add(e1);
        Polynomial a = v.mult(aKey).add(e2);

        if (debug) {
            System.out.println("final b");
            b.debugPrint();
            System.out.println("");
            System.out.println("final a");
            a.debugPrint();
        }

        res.init(b, a, src.getLevel(), src.getScale());
    }
}
