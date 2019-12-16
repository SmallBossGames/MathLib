package smallBossMathLib.shared

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction

class UnivariateDifferentiableKotlinFunction(
    val diffFunc: (t: DerivativeStructure) -> DerivativeStructure,
    val func: (x: Double) -> Double)
    : UnivariateDifferentiableFunction {

    override fun value(t: DerivativeStructure?): DerivativeStructure {
        requireNotNull(t)
        return diffFunc(t)
    }

    override fun value(x: Double): Double {
        return func(x)
    }
}