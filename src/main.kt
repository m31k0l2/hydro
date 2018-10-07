import kotlin.math.*

// n - обводненность
class Hydro(val ro1: Double, val ro2: Double, mu1: Double, mu2: Double, val Gp: Double, val n: Double, t: Double) {
    private val To = 273.15
    private val Pa = 101000.0
    private val T = To + t
    private val g = 9.81
    private val muOt = mu2/mu1
    private val v1 = mu1/ro1
    private val v2 = mu2/ro2
    private val Ke = 0.2 // https://studfiles.net/preview/2474060/page:5/
    private val sigma = 30/1000.0
    val A = 1.3342
    val B = -0.3342
    val C = 0.27
    val zg = 0.8
    fun Q2(Gsv: Double, Qn: Double, P: Double) = Gsv*Qn*T*Pa*zg/(To*P)
    fun Gsv(P: Double) = (A+B*Math.pow(P/Pa, C))*Gp
    fun dPPerDl(lambdaSm: Double, wSm: Double, roPs: Double, roFi: Double, sinA: Double, D: Double) = lambdaSm*wSm*wSm*roPs/(2*D) + roFi*g*sinA
    fun roPs(b1: Double, fi1: Double, b2: Double, fi2: Double) = b1*b1*ro1/fi1 + b2*b2*ro2/fi2
    fun roFi(fi1: Double, fi2: Double) = fi1*ro1 + fi2*ro2
    fun vStar(Re2: Double, Frsm: Double, b1: Double) =
            10*(0.82-0.0017*Math.pow(muOt, 0.6))*Math.pow((Re2*Frsm*ro2/(ro1-ro2)), -1.0/3)*Math.exp((8+62*muOt)*b1)
    fun Re2(wSm: Double, D: Double) = wSm*D/v2
    fun FrSm(wSm: Double, D: Double) = wSm*wSm/(g*D)
    fun FrStar(sinA: Double, lambda2: Double, b1: Double, b2: Double) = (0.2+ 2*sinA/lambda2)* exp(-2.5*b2)/b1/b1
    fun lambda2(sinA: Double, D: Double, w: Double) = 2*sinA*g*D/w
    fun fi1kv(fi1Star: Double, b1: Double, A: Double) = fi1Star/(1+200*b1) + 55*Math.pow(b1, 0.5)/A
    fun A(Re1: Double, Frsm: Double) = Math.pow((Re1*Frsm*ro2/(ro1-ro2)), 1/3.0)
    fun fi1Star(Wstar: Double, A: Double) = if (Math.round(Wstar*10)/10.0 >= 3.3) 0.0 else 0.0053*(3.3-Wstar)/A
    fun Wstar(wSm: Double, sinA: Double) = Math.pow(wSm*((ro1-ro2)/(sigma*g*sinA)), 0.25)*Math.pow(ro2/ro1, 0.5)
    fun fi1kn(A: Double, sinA: Double, b1: Double) = 55*(1-1/(1+3.84*Math.pow(10.0, -6.0)*A*A*A*Math.pow(Math.abs(sinA), -1.66)))*Math.pow(b1, 0.5)/A
    fun fi1pv(b2: Double, Frsm: Double) = 1 - K1*b2*(1- exp(-4.4*Math.pow(Frsm/Fra, 0.5)))
    val K1 = if (round(muOt*100)/100.0 <= 0.01) 0.35 + 1.4*Math.pow(muOt, 0.25) else 0.77 + 0.23*Math.pow(muOt, 0.5)
    val Fra = if (round(muOt*1000)/1000.0 >= 0.001) 9.8*Math.pow(muOt, 0.1) else 1150*Math.pow(muOt, 0.79)
    fun fi1pn(b2: Double) = 1 - K1*b2
    fun lambdaSm(lambda0: Double, psi: Double) = lambda0*psi
    fun lambda0(Resm: Double, D: Double) = 0.067*(158/Resm + 2*Ke/D)
    fun vsm(b1: Double, b2: Double, isRound: Boolean) = if (isRound) v1 else 1/(b1/v1 + b2/v2)
    fun psi(Re2: Double, Frsm: Double, b1: Double, b2: Double, isRound: Boolean) = if (isRound)
        1 + 0.031*Math.pow((Re2*Frsm*(ro1-ro2)/ro2), 1/3.0)* exp(-15*(ro2/ro1+b1))*Math.pow(b1, 0.5) else
        (1-0.78*b2*(1- exp(-4.4*Frsm/Fra)) - 0.22*b2*(1- exp(-15*ro2/ro2)))/b1
    fun Resm(wSm: Double, D: Double, vsm: Double) = wSm*D/vsm
    fun wSm(Qsm: Double, D: Double) = Qsm/Ftr(D)
    fun Ftr(D: Double) = Math.PI*D*D/4
    fun b2(Q1: Double, Q2: Double) = Q2/(Q1+Q2)
    fun isUp(H1: Double, H2: Double) = H1 < H2
    fun Re1(wSm: Double, D: Double) = wSm*D/v1
    var info = ""

    fun calculate(D: Double, H1: Double, H2: Double, L: Double, Qsm: Double, Pk: Double): Double {
        info += "Q = $Qsm\n"
        info += "D = $D\n"
        info += "Pk = $Pk\n"
        val Qn = Qsm*(1-n)
        info += "Qn = $Qn\n"
        val Gsv = Gsv(Pk)
        info += "Gsv = $Gsv\n"
        val Q2 = Q2(Gsv, Qn, Pk)
        info += "Q2 = $Q2\n"
        val Q1 = (Qsm*ro1 - Q2*ro2)/ro1
        info += "Q1 = $Q1\n"
        val wSm = wSm(Qsm, D)
        info += "wSm = $wSm\n"
        val b2 = b2(Q1, Q2)
        info += "b2 = $b2\n"
        val b1 = 1 - b2
        info += "b1 = $b1\n"
        val Re2 = Re2(wSm, D)
        info += "Re2 = $Re2\n"
        val Frsm = FrSm(wSm, D)
        info += "Frsm = $Frsm\n"
        val vStar = vStar(Re2, Frsm, b1)
        info += "V* = $vStar\n"
        val sinA = (H2-H1)/L
        info += "sinA = $sinA\n"
        val rezhim = if (Math.round(vStar*10)/10.0 <= 1.0) 0 else {
            val lambda2 = lambda2(sinA, D, wSm)
            info += "lambda2 = $lambda2\n"
            val FrStar = FrStar(sinA, lambda2, b1, b2)
            info += "Fr* = $FrStar\n"
            if (Frsm >= FrStar) 1 else 2
        }
        info += (if (rezhim == 0) "Кольцевой" else if (rezhim == 1) "Пробковый" else "Расслоенный") + "\n"
        val vsm = vsm(b1, b2, rezhim == 0)
        info += "vsm = $vsm\n"
        val Resm = Resm(wSm, D, vsm)
        info += "Resm = $Resm\n"
        val psi = psi(Re2, Frsm, b1, b2, rezhim == 0)
        info += "psi = $psi\n"
        if (rezhim < 2) {
            val lambda0 = lambda0(Resm, D)
            info += "lambda0 = $lambda0\n"
            val lambdaSm = lambdaSm(lambda0, psi)
            info += "lambdaSm = $lambdaSm\n"
            val Re1 = Re1(wSm, D)
            info += "Re1 = $Re1\n"
            val A = A(Re1, Frsm)
            info += "A = $A\n"
            val fi1 = if (rezhim == 0) {
                if (isUp(H1, H2)) {
                    val Wstar = Wstar(wSm, sinA)
                    info += "W* = $Wstar\n"
                    val fi1Star = fi1Star(Wstar, A)
                    info += "fi1* = $fi1Star\n"
                    fi1kv(fi1Star, b1, A)
                } else {
                    fi1kn(A, sinA, b1)
                }
            } else {
                if (isUp(H1, H2)) {
                    fi1pv(b2, Frsm)
                } else {
                    fi1pn(b2)
                }
            }
            info += "fi1 = $fi1\n"
            val fi2 = 1 - fi1
            info += "fi2 = $fi2\n"
            val roPs = roPs(b1, fi1, b2, fi2)
            info += "roPs = $roPs\n"
            val roFi = roFi(fi1, fi2)
            info += "roFi = $roFi\n"
            val dpPerDl = dPPerDl(lambdaSm, wSm, roPs, roFi, sinA, D)
            info += "dp/dL = $dpPerDl\n"
            val dp = dpPerDl * L
            info += "dP = $dp\n"
            info += "Pn = ${Pk + dp}"
            return Pk + dp
        } else {
            val fi2 = 2.63*Math.pow((Frsm*ro2/(ro1-ro2)*b1*b1/(0.02+Math.pow(sinA, 0.5))), 1/3.0)
            info += "fi2 = $fi2\n"
            var tetta = 0.0
            var fi2Star = 0.0
            var min = Double.MAX_VALUE
            for (i in 0..180) {
                val x = i/Math.PI
                val fi = tetta - sin(tetta)*cos(tetta)
                if (abs(fi - fi2) < min) {
                    min = abs(fi2 - fi)
                    tetta = x
                    fi2Star = fi
                }
            }
            info += "fi2* = $fi2Star\n"
            info += "tetta = $tetta\n"
            val Dr = fi2*Math.PI*D/tetta
            info += "Dr = $Dr\n"
            val dpPerDl = lambda2(sinA, D, wSm)*ro2*wSm*wSm/2/Dr+ro2*g*sinA
            info += "dpPerDl = $dpPerDl\n"
            val dp = dpPerDl * L
            info += "dP = $dp\n"
            info += "Pn = ${Pk + dp}"
            return Pk + dp
        }
    }
}

fun main(args: Array<String>) {
    val h = Hydro(1090.9, 1.7, 0.0044, 0.00001, 22.0, 0.7, 20.0)
    val p = h.calculate(0.3, 205.0, 200.0, 500.0, 13634.9/24/3600, 0.3*Math.pow(10.0, 6.0))
    println(p/Math.pow(10.0, 6.0))
}