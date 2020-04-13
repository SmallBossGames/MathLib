package smallBossMathLib.implicitDifferentialEquations.exceptions

class ExceedingLimitEvaluationsException : Exception() {
    override val message: String?
        get() = "Maximum evaluations exceeded"
}