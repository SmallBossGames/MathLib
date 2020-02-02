package smallBossMathLib.explicitDifferentialEquations

abstract class RungeKuttaIntegratorBase {
    private val stepHandlers = ArrayList<(t: Double, y: DoubleArray) -> Unit>()

    fun addStepHandler(handler: (t: Double, y: DoubleArray) -> Unit) {
        stepHandlers += handler
    }

    fun removeStepHandler(handler: (t: Double, y: DoubleArray) -> Unit) {
        stepHandlers.remove(handler)
    }

    fun clearStepHandlers() {
        stepHandlers.clear()
    }

    protected fun executeStepHandlers(t: Double, y: DoubleArray) {
        for (handler in stepHandlers)
            handler(t, y)
    }
}