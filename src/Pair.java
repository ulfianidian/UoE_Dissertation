public class Pair implements Comparable<Pair>{
    public final int index;
    public final double value;

    public Pair(int index, double value){
        this.index = index;
        this.value = value;
    }

    @Override
    public int compareTo(Pair other){
        return -1 * Double.valueOf(this.value).compareTo(other.value);
    }

    public int getKey() {
        return this.index;
    }

    public double getValue() { return this.value; }
}
