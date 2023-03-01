package ckks;

import java.util.Arrays;

public class Evaluator {
    final boolean debug = false;

    private Context context;

    private PublicKeys publicKey;

    public Evaluator(Context context, PublicKeys publicKey) {
        this.context = context;
        this.publicKey = publicKey;
    }

    // add
    public Plaintext add(Plaintext a, Plaintext b) {
        Plaintext res = new Plaintext(a);
        add_inplace(res, b);
        return res;
    }

    public void add_inplace(Plaintext dest, Plaintext other) {
        assertNotEmpty(dest);
        assertNotEmpty(other);

        Polynomial mRes = dest.getM();
        Polynomial mOther = other.getM();

        mRes.add_inplace(mOther);
    }

    public Ciphertext add(Ciphertext a, Plaintext b) {
        Ciphertext res = new Ciphertext(a);
        add_inplace(res, b);
        return res;
    }

    public void add_inplace(Ciphertext dest, Plaintext other) {
        assertNotEmpty(dest);
        assertNotEmpty(other);

        if (dest.isLinear()) {
            Polynomial bRes = dest.getB();
            Polynomial m = other.getM();

            bRes.add_inplace(m);
        } else {
            throw new UnsupportedOperationException("not implemented");
        }
    }

    public Ciphertext add(Ciphertext a, Ciphertext b) {
        Ciphertext res = new Ciphertext(a);
        add_inplace(res, b);
        return res;
    }

    public void add_inplace(Ciphertext dest, Ciphertext other) {
        assertNotEmpty(dest);
        assertNotEmpty(other);
        assertSameRelinearizationStatus(dest, other);

        if (dest.isLinear()) {
            Polynomial aRes = dest.getA();
            Polynomial bRes = dest.getB();
            Polynomial aOther = other.getA();
            Polynomial bOther = other.getB();

            aRes.add_inplace(aOther);
            bRes.add_inplace(bOther);
        } else {
            Polynomial d0Res = dest.getD0();
            Polynomial d1Res = dest.getD1();
            Polynomial d2Res = dest.getD2();
            Polynomial d0Other = other.getD0();
            Polynomial d1Other = other.getD1();
            Polynomial d2Other = other.getD2();

            d0Res.add_inplace(d0Other);
            d1Res.add_inplace(d1Other);
            d2Res.add_inplace(d2Other);
        }
    }

    // mult
    public Plaintext mult(Plaintext p, int integer) {
        Plaintext res = new Plaintext(p);
        mult_inplace(res, integer);
        return res;
    }

    public void mult_inplace(Plaintext dest, int integer) {
        assertNotEmpty(dest);

        Polynomial mRes = dest.getM();

        mRes.mult_inplace(integer);
    }

    public Plaintext mult(Plaintext a, Plaintext b) {
        Plaintext res = new Plaintext(a);
        mult_inplace(res, b);
        return res;
    }

    public void mult_inplace(Plaintext dest, Plaintext other) {
        assertNotEmpty(dest);
        assertNotEmpty(other);

        assertNonZeroLevel(dest);
        assertNonZeroLevel(other);

        Polynomial mRes = dest.getM();
        Polynomial mOther = other.getM();

        mRes.mult_inplace(mOther);

        dest.setScale(dest.getScale() * other.getScale());
    }

    public Ciphertext mult(Ciphertext c, int integer) {
        Ciphertext res = new Ciphertext(c);
        mult_inplace(res, integer);
        return res;
    }

    public void mult_inplace(Ciphertext dest, int integer) {
        assertNotEmpty(dest);

        Polynomial aRes = dest.getA();
        Polynomial bRes = dest.getB();

        aRes.mult_inplace(integer);
        bRes.mult_inplace(integer);
    }

    public Ciphertext mult(Ciphertext a, Plaintext b) {
        Ciphertext res = new Ciphertext(a);
        mult_inplace(res, b);
        return res;
    }

    public void mult_inplace(Ciphertext dest, Plaintext other) {
        assertNotEmpty(dest);
        assertNotEmpty(other);

        assertNonZeroLevel(dest);
        assertNonZeroLevel(other);

        assertRelinearized(dest);

        Polynomial aRes = dest.getA();
        Polynomial bRes = dest.getB();
        Polynomial m = other.getM();

        bRes.mult_inplace(m);
        aRes.mult_inplace(m);

        dest.setScale(dest.getScale() * other.getScale());
    }

    public Ciphertext mult(Ciphertext a, Ciphertext b) {
        Ciphertext res = new Ciphertext(a);
        mult_inplace(res, b);
        return res;
    }

    public void mult_inplace(Ciphertext dest, Ciphertext other) {
        assertNotEmpty(dest);
        assertNotEmpty(other);

        assertNonZeroLevel(dest);
        assertNonZeroLevel(other);

        assertRelinearized(dest);
        assertRelinearized(other);

        Polynomial aRes = dest.getA();
        Polynomial bRes = dest.getB();
        Polynomial aOther = other.getA();
        Polynomial bOther = other.getB();

        Polynomial d0 = bRes.mult(bOther);
        Polynomial d1 = bRes.mult(aOther).add(aRes.mult(bOther));
        Polynomial d2 = aRes.mult(aOther);

        dest.afterMult(d0, d1, d2, other.getScale() * dest.getScale());
    }

    // relinearize
    public Ciphertext relinearize(Ciphertext src) {
        Ciphertext res = new Ciphertext(src);
        relinearize_inplace(res);
        return res;
    }

    public void relinearize_inplace(Ciphertext src) {
        // already linear
        if (src.getD2() == null)
            return;

        long[][] d2 = src.getD2().getCrt();

        if (debug) {
            System.out.println("Evaluator.relinearize_inplace");
            System.out.println("src D2");
            System.out.println(Arrays.deepToString(d2));
            System.out.println();
        }

        Maths.nttInverse_inplace(context, d2, src.getLevel());

        if (debug) {
            System.out.println("D2 after NTT inverse");
            System.out.println(Arrays.deepToString(d2));
            System.out.println();
        }

        long[][] d2ModedUp = Maths.modUp(context, d2, src.getLevel());

        if (debug) {
            System.out.println("D2 after mod up");
            System.out.println(Arrays.deepToString(d2ModedUp));
            System.out.println();
        }

        long[][] c0Coeffs = new long[d2ModedUp.length][d2ModedUp[0].length];
        long[][] c1Coeffs = new long[d2ModedUp.length][d2ModedUp[0].length];

        long[][] relinKeyB = publicKey.getRelinKeyB();
        long[][] relinKeyA = publicKey.getRelinKeyA();

        if (debug) {
            System.out.println("relinKeyB");
            System.out.println(Arrays.deepToString(relinKeyB));
            System.out.println();
            System.out.println("relinKeyA");
            System.out.println(Arrays.deepToString(relinKeyA));
            System.out.println();
        }

        // mult by relinearization key
        for (int primeIdx = 0; primeIdx <= src.getLevel(); primeIdx++)
            c0Coeffs[primeIdx] = Maths.multPolynominalsMod(d2ModedUp[primeIdx], relinKeyB[primeIdx],
                    context.primes[primeIdx]);
        for (int tempPrimeIdx = 0; tempPrimeIdx < context.tempPrimes.length; tempPrimeIdx++)
            c0Coeffs[context.primes.length + tempPrimeIdx] = Maths.multPolynominalsMod(
                    d2ModedUp[context.primes.length + tempPrimeIdx], relinKeyB[context.primes.length + tempPrimeIdx],
                    context.tempPrimes[tempPrimeIdx]);

        for (int primeIdx = 0; primeIdx <= src.getLevel(); primeIdx++)
            c1Coeffs[primeIdx] = Maths.multPolynominalsMod(d2ModedUp[primeIdx], relinKeyA[primeIdx],
                    context.primes[primeIdx]);
        for (int tempPrimeIdx = 0; tempPrimeIdx < context.tempPrimes.length; tempPrimeIdx++)
            c1Coeffs[context.primes.length + tempPrimeIdx] = Maths.multPolynominalsMod(
                    d2ModedUp[context.primes.length + tempPrimeIdx], relinKeyA[context.primes.length + tempPrimeIdx],
                    context.tempPrimes[tempPrimeIdx]);

        if (debug) {
            System.out.println("c0Coeffs");
            System.out.println(Arrays.deepToString(c0Coeffs));
            System.out.println();
            System.out.println("c1Coeffs");
            System.out.println(Arrays.deepToString(c1Coeffs));
            System.out.println();
        }

        long[][] c0ntt = Maths.ntt(context, Maths.modDown(context, c0Coeffs, src.getLevel()), src.getLevel());
        long[][] c1ntt = Maths.ntt(context, Maths.modDown(context, c1Coeffs, src.getLevel()), src.getLevel());

        if (debug) {
            System.out.println("C0 after mod down and NTT");
            System.out.println(Arrays.deepToString(c0ntt));
            System.out.println();
            System.out.println("C1 after mod down and NTT");
            System.out.println(Arrays.deepToString(c1ntt));
            System.out.println();
        }

        Polynomial c0 = new Polynomial(context, c0ntt);
        Polynomial c1 = new Polynomial(context, c1ntt);

        // add c0 to d0 and c1 to d1
        Polynomial b = src.getD0().add(c0);
        Polynomial a = src.getD1().add(c1);

        if (debug) {
            System.out.println("final b");
            b.debugPrint();
            System.out.println();
            System.out.println("final a");
            a.debugPrint();
            System.out.println();
        }

        src.afterRelinearization(b, a);
    }

    // rescale
    public Ciphertext rescale(Ciphertext src) {
        Ciphertext res = new Ciphertext(src);
        rescale_inplace(res);
        return res;
    }

    public void rescale_inplace(Ciphertext src) {
        assertRelinearized(src);
        assertNonZeroLevel(src);

        // already scaled correctly
        if (src.getScale() < 2 * context.defaultScale)
            return;

        if (debug) {
            System.out.println("Evaluator.rescale_inplace");
            System.out.println("src B");
            System.out.println(Arrays.deepToString(src.getB().getCrt()));
            System.out.println("src A");
            System.out.println(Arrays.deepToString(src.getA().getCrt()));
            System.out.println("src level= " + src.getLevel() + ", src scale= " + src.getScale());
            System.out.println();
        }

        long[][] c0coeffs = Maths.nttInverse(context, src.getB().getCrt(),
                src.getLevel());
        long[][] c1coeffs = Maths.nttInverse(context, src.getA().getCrt(),
                src.getLevel());

        if (debug) {
            System.out.println("c0 coeffs");
            System.out.println(Arrays.deepToString(c0coeffs));
            System.out.println("c1 coeffs");
            System.out.println(Arrays.deepToString(c1coeffs));
            System.out.println("Top level prime= " + context.primes[src.getLevel()]);
            System.out.println();
        }

        for (int primeIdx = 0; primeIdx < src.getLevel(); primeIdx++) {
            if (debug) {
                System.out.println("Current prime= " + context.primes[primeIdx]
                        + ", top level prime inverse modulo current prime= "
                        + context.otherPrimesInv[primeIdx][src.getLevel()] + "\n");
            }

            for (int coeffIdx = 0; coeffIdx < context.slots * 2; coeffIdx++) {
                c0coeffs[primeIdx][coeffIdx] = Maths.mod(
                        c0coeffs[primeIdx][coeffIdx] - c0coeffs[src.getLevel()][coeffIdx], context.primes[primeIdx]);
                c0coeffs[primeIdx][coeffIdx] = Maths.modMult(c0coeffs[primeIdx][coeffIdx],
                        context.otherPrimesInv[primeIdx][src.getLevel()], context.primes[primeIdx]);

                c1coeffs[primeIdx][coeffIdx] = Maths.mod(
                        c1coeffs[primeIdx][coeffIdx] - c1coeffs[src.getLevel()][coeffIdx], context.primes[primeIdx]);
                c1coeffs[primeIdx][coeffIdx] = Maths.modMult(c1coeffs[primeIdx][coeffIdx],
                        context.otherPrimesInv[primeIdx][src.getLevel()], context.primes[primeIdx]);
            }
        }

        if (debug) {
            System.out.println("c0 coeffs rescaled");
            System.out.println(Arrays.deepToString(c0coeffs));
            System.out.println("c1 coeffs rescaled");
            System.out.println(Arrays.deepToString(c1coeffs));
            System.out.println();
        }

        long[][] c0ntt = Maths.ntt(context, c0coeffs, src.getLevel() - 1);
        long[][] c1ntt = Maths.ntt(context, c1coeffs, src.getLevel() - 1);

        Polynomial b = new Polynomial(context, c0ntt);
        Polynomial a = new Polynomial(context, c1ntt);

        if (debug) {
            System.out.println("Final b");
            b.debugPrint();
            System.out.println("Final a");
            a.debugPrint();
            System.out.println("New level= " + (src.getLevel() - 1) + ", new scale= "
                    + (src.getScale() / context.primes[src.getLevel()]));
            System.out.println();
        }

        src.afterRescale(b, a, src.getScale() / context.primes[src.getLevel()], src.getLevel() - 1);
    }

    private void assertNotEmpty(Ciphertext c) {
        if (c.isEmpty())
            throw new IllegalStateException("Ciphertext is empty.");
    }

    private void assertNotEmpty(Plaintext p) {
        if (p.isEmpty())
            throw new IllegalStateException("Plaintext is empty.");
    }

    private void assertNonZeroLevel(Ciphertext c) {
        if (c.getLevel() == 0)
            throw new IllegalStateException("Ciphertext is on level 0.");
    }

    private void assertNonZeroLevel(Plaintext p) {
        if (p.getLevel() == 0)
            throw new IllegalStateException("Plaintext is on level 0.");
    }

    private void assertRelinearized(Ciphertext c) {
        if (c.isEmpty())
            throw new IllegalStateException("Ciphertext is not relinearized.");
    }

    private void assertSameRelinearizationStatus(Ciphertext a, Ciphertext b) {
        if (a.isLinear() != b.isLinear())
            throw new IllegalStateException("Both ciphertexts must be either relinearized or not relinearized");
    }
}
