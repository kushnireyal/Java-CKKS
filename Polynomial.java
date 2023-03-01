package ckks;

import java.util.Arrays;

public class Polynomial {
    final static boolean debug = false;

    private Context context;

    private long[][] crt;

    // c'tors
    public Polynomial() {
    }

    public Polynomial(Context context, long[][] crt) {
        this.context = context;
        this.crt = crt;
    }

    public Polynomial(Polynomial poly) {
        this.context = poly.context;
        this.crt = new long[context.primes.length][context.slots * 2];
        for (int i = 0; i < context.primes.length; i++)
            for (int j = 0; j < context.slots * 2; j++)
                this.crt[i][j] = poly.crt[i][j];
    }

    // methods
    public Polynomial add(Polynomial that) {
        Polynomial res = new Polynomial(this);
        res.add_inplace(that);
        return res;
    }

    public void add_inplace(Polynomial that) {
        for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
            for (int i = 0; i < context.slots * 2; i++)
                crt[primeIdx][i] = Maths.mod(crt[primeIdx][i] + that.crt[primeIdx][i], context.primes[primeIdx]);
    }

    public Polynomial sub(Polynomial that) {
        Polynomial res = new Polynomial(this);
        res.sub_inplace(that);
        return res;
    }

    public void sub_inplace(Polynomial that) {
        for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
            for (int i = 0; i < context.slots * 2; i++)
                crt[primeIdx][i] = Maths.mod(crt[primeIdx][i] - that.crt[primeIdx][i], context.primes[primeIdx]);
    }

    public Polynomial mult(Polynomial that) {
        Polynomial res = new Polynomial(this);
        res.mult_inplace(that);
        return res;
    }

    public void mult_inplace(Polynomial that) {
        if (debug) {
            System.out.println("Polynomial mult_inplace:");
            System.out.println("this crt:");
            System.out.println(Arrays.deepToString(crt));
            System.out.println("that crt:");
            System.out.println(Arrays.deepToString(that.crt));
        }

        for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
            for (int i = 0; i < context.slots * 2; i++) {
                if (debug)
                    System.out.println("calculating " + crt[primeIdx][i] + " * " + that.crt[primeIdx][i] + " mod "
                            + context.primes[primeIdx]);
                crt[primeIdx][i] = Maths.modMult(crt[primeIdx][i], that.crt[primeIdx][i], context.primes[primeIdx]);
                if (debug)
                    System.out.println("result= " + crt[primeIdx][i]);
            }

        if (debug) {
            System.out.println("res crt:");
            System.out.println(Arrays.deepToString(crt));
            System.out.println();
        }
    }

    public Polynomial mult(int scalar) {
        Polynomial res = new Polynomial(this);
        res.mult_inplace(scalar);
        return res;
    }

    public void mult_inplace(int scalar) {
        for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
            for (int i = 0; i < context.slots * 2; i++)
                crt[primeIdx][i] = Maths.modMult(crt[primeIdx][i], scalar, context.primes[primeIdx]);
    }

    public Polynomial square() {
        Polynomial res = new Polynomial(this);
        res.square_inplace();
        return res;
    }

    public void square_inplace() {
        mult_inplace(this);
    }

    public long[][] getCrt() {
        return crt;
    }

    public long[] serialize() {
        long[] res = new long[2 + crt.length * crt[0].length];

        int idx = 0;
        res[idx++] = crt.length;
        res[idx++] = crt[0].length;

        for (int i = 0; i < crt.length; i++)
            for (int j = 0; j < crt[0].length; j++)
                res[idx++] = crt[i][j];

        return res;
    }

    public void deserialize(Context context, long[] serialization) {
        this.context = context;

        int idx = 0;

        crt = new long[(int) serialization[idx++]][(int) serialization[idx++]];

        for (int i = 0; i < crt.length; i++)
            for (int j = 0; j < crt[0].length; j++)
                crt[i][j] = serialization[idx++];
    }

    public void debugPrint() {
        System.out.println(Arrays.deepToString(crt));
    }
}
