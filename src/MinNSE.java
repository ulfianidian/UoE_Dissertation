public class MinNSE {
    private double mValue;
    private double yValue;
    private double leftAllot;

    public MinNSE(double mValue, double yValue, double leftAllot){
        this.mValue = mValue;
        this.yValue = yValue;
        this.leftAllot = leftAllot;
    }

    public double getmValue() { return mValue; }

    public double getyValue() { return yValue; }

    public double getLeftAllot() { return leftAllot; }

    public void setmValue(double newMValue) { this.mValue = newMValue; }

    public void setyValue(double newYValue) { this.yValue = newYValue; }

    public void setLeftAllot(double newLeftAllot) { this.leftAllot = newLeftAllot; }
}
