package smallBossMathLib.shared

abstract class ExplicitIntegrator : Integrator() {
    private val stepHandlers = HashSet<(time:Double, y:DoubleArray, info:IExplicitMethodStepInfo) -> Unit>()

    protected fun executeStepHandlers(time:Double, y:DoubleArray, info:IExplicitMethodStepInfo) {
        for (handler in stepHandlers)
            handler(time, y, info)
    }

    fun addStepHandler(handler: (time:Double, y:DoubleArray, info:IExplicitMethodStepInfo) -> Unit) {
        stepHandlers += handler
    }

    fun removeStepHandler(handler: (time:Double, y:DoubleArray, info:IExplicitMethodStepInfo) -> Unit) {
        stepHandlers -= handler
    }
}