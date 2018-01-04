import gngs.Pedigree
import gngs.Variant

/**
 * Represents a constraint in a genetic model
 * 
 * @author simon.sadedin@mcri.edu.au
 */
interface Constraint {
    
    boolean test(Variant v, List<Pedigree> pedigrees);

}
