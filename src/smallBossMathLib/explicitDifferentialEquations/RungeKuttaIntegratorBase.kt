package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.differentialEquations.IFirstOrderIntegrator

abstract class RungeKuttaIntegratorBase {
    var isEvaluationControlEnabled: Boolean = false
        private set

    var maxEvaluationCount: Int = 0
        private set

    private val stepHandlers = HashSet<(t: Double, y: DoubleArray) -> Unit>()

    fun addStepHandler(handler: (t: Double, y: DoubleArray) -> Unit) {
        stepHandlers += handler
    }

    fun removeStepHandler(handler: (t: Double, y: DoubleArray) -> Unit) {
        stepHandlers -= handler
    }

    fun clearStepHandlers() {
        stepHandlers.clear()
    }

    protected fun executeStepHandlers(t: Double, y: DoubleArray) {
        for (handler in stepHandlers)
            handler(t, y)
    }

    //Evaluation control
    fun enableEvaluationCountCheck(maxEvaluations: Int){
        isEvaluationControlEnabled = true
        maxEvaluationCount = maxEvaluations
    }

    fun disableEvaluationCountCheck(){
        isEvaluationControlEnabled = false
        maxEvaluationCount = 0
    }

    protected fun isNextEvaluationAllow(evaluations: Int) : Boolean =
        !isEvaluationControlEnabled || evaluations < maxEvaluationCount
}