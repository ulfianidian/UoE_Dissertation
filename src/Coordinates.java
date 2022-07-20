public class Coordinates {
    private final Integer column, row;

    public Coordinates(int column, int row){
        this.column = column;
        this.row = row;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Coordinates))
            return false;
        Coordinates other = (Coordinates) obj;
        return column == other.column && row == other.row;
    }

    @Override
    public int hashCode() {
        return 31 * column.hashCode() + row.hashCode();
    }
}

//public class SparseMatrix{
//
//    private Map<Coordinates, MinNSE> minNSEMap = new HashMap<>();
//    private int maxColumns, maxRows;
//
//    public SparseMatrix(int maxColumns, int maxRows){
//        this.maxColumns = maxColumns;
//        this.maxRows = maxRows;
//    }
//
//    public void SetValue(int column, int row, MinNSE minNSE){
//        if(column > maxColumns || row > maxRows)
//            throw new RuntimeException("Position out of bounds");
//        minNSEMap.put(new Coordinates(column, row), minNSE);
//    }
//
//    public Map<Coordinates, MinNSE> getMapping() { return minNSEMap; }
//}
