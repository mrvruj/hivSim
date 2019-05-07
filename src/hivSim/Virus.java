package hivSim;

/**
 * An instance is used as the founder (root) and descendants (nodes) of a
 * VirusTree object.
 */
public class Virus<T> {
    private T data;
    private double burstday;

    /** Constructor: an instance of Virus with specified data and burstday.
     *
     * @param data the information (of any type) stored in the Virus object.
     * @param burstday a double which represents the time the Virus
     *        object and it's corresponding VirusTree release new
     *        Viruses/VirusTrees. Once burst, Viruses should not generate
     *        new Viruses/VirusTrees. Must be >= 0.0.
     */
    public Virus(T data, double burstday) {
        if (burstday < 0){
            throw new IllegalArgumentException("Burstday must be >= 0.0.");
        }
        this.data = data;
        this.burstday = burstday;
    }

    public void setData(T data) {
        this.data = data;
    }

    /** @return the data stored in this Virus. */
    public T getData() {
        return this.data;
    }

    public void setBurstday(double burstday) {
        if (burstday < 0){
            throw new IllegalArgumentException("Burstday must be >= 0.0.");
        }
        this.burstday = burstday;
    }

    /** @return the burstday of this Virus. */
    public double getBurstday() {
        return this.burstday;
    }

    @Override
    public String toString() {
        if (data == null){return "";}
        return this.data.toString();
    }
}
