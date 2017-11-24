package com.pescamillam.fem.element;

public class Cst {
    private Point i;
    private Point j;
    private Point m;
    public Cst(Point i, Point j, Point m) {
        super();
        this.i = i;
        this.j = j;
        this.m = m;
    }
    public Point getI() {
        return i;
    }
    public void setI(Point i) {
        this.i = i;
    }
    public Point getJ() {
        return j;
    }
    public void setJ(Point j) {
        this.j = j;
    }
    public Point getM() {
        return m;
    }
    public void setM(Point m) {
        this.m = m;
    }
}
