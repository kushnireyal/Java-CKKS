package ckks;

import java.math.BigInteger;

public class Complex {

    public double re, im;

    // c'tors
    public Complex(double re, double im) {
        this.re = re;
        this.im = im;
    }

    public Complex(double re) {
        this(re, 0);
    }

    public Complex(BigInteger re) {
        this(re.doubleValue(), 0);
    }

    public Complex() {
        this(0, 0);
    }

    public Complex(Complex src) {
        this(src.re, src.im);
    }

    // basic arithmetics
    public Complex add(Complex that) {
        return new Complex(this.re + that.re, this.im + that.im);
    }

    public void add_inplace(Complex that) {
        this.re += that.re;
        this.im += that.im;
    }

    public Complex sub(Complex that) {
        return new Complex(this.re - that.re, this.im - that.im);
    }

    public void sub_inplace(Complex that) {
        this.re -= that.re;
        this.im -= that.im;
    }

    public Complex mult(Complex that) {
        return new Complex(this.re * that.re - this.im * that.im, this.re * that.im + this.im * that.re);
    }

    public void mult_inplace(Complex that) {
        double new_re = this.re * that.re - this.im * that.im;
        double new_im = this.re * that.im + this.im * that.re;

        this.re = new_re;
        this.im = new_im;
    }

    public Complex divide(Complex that) {
        Complex res = new Complex(this.mult(that.conj()));

        double thatNormSquared = that.mult(that.conj()).re;

        res.re /= thatNormSquared;
        res.im /= thatNormSquared;

        return res;
    }

    public void divide_inplace(Complex that) {
        double thatNormSquared = that.mult(that.conj()).re;

        this.mult_inplace(that.conj());

        this.re /= thatNormSquared;
        this.im /= thatNormSquared;
    }

    public Complex conj() {
        return new Complex(this.re, -this.im);
    }

    public void conj_inplace() {
        this.im = -this.im;
    }

    public double norm() {
        return Math.sqrt(this.re * this.re + this.im * this.im);
    }

    public boolean equal(Complex other) {
        return this.re == other.re && this.im == other.im;
    }

    // print
    @Override
    public String toString() {
        return "(" + this.re + ", " + this.im + ")";
    }

    // vector operations
    static public Complex hermitianProduct(Complex[] a, Complex[] b) {
        if (a.length != b.length)
            throw new IllegalArgumentException("Both vectors must be of the same size.");

        Complex res = new Complex(0, 0);

        for (int i = 0; i < a.length; i++)
            res.add_inplace(a[i].mult(b[i].conj()));

        return res;
    }

    static public Complex[] multVectorByScalar(Complex[] vec, Complex scalar) {
        Complex[] res = deepClone(vec);

        for (int i = 0; i < res.length; i++)
            res[i].mult_inplace(scalar);

        return res;
    }

    static public Complex[] multVectorByScalar(Complex[] vec, double scalar) {
        return multVectorByScalar(vec, new Complex(scalar));
    }

    static public void multVectorByScalar_inplace(Complex[] vec, Complex scalar) {
        for (int i = 0; i < vec.length; i++)
            vec[i].mult_inplace(scalar);
    }

    static public void multVectorByScalar_inplace(Complex[] vec, double scalar) {
        multVectorByScalar_inplace(vec, new Complex(scalar));
    }

    static public Complex[] addVectors(Complex[] a, Complex[] b) {
        if (a.length != b.length)
            throw new IllegalArgumentException("Both vectors must be of the same size.");

        Complex[] res = deepClone(a);

        for (int i = 0; i < a.length; i++)
            res[i].add_inplace(b[i]);

        return res;
    }

    static public void addVectors_inplace(Complex[] dest, Complex[] other) {
        if (dest.length != other.length)
            throw new IllegalArgumentException("Both vectors must be of the same size.");

        for (int i = 0; i < dest.length; i++)
            dest[i].add_inplace(other[i]);
    }

    static public Complex[] deepClone(Complex[] src) {
        Complex[] clone = new Complex[src.length];

        for (int i = 0; i < src.length; i++)
            clone[i] = new Complex(src[i]);

        return clone;
    }
}
