import javafx.fxml.FXML
import javafx.scene.control.TabPane
import javafx.scene.control.TextArea
import javafx.scene.control.TextField

class Form {
    @FXML lateinit var D: TextField
    @FXML lateinit var H1: TextField
    @FXML lateinit var H2: TextField
    @FXML lateinit var L: TextField
    @FXML lateinit var Q: TextField
    @FXML lateinit var Pk: TextField
    @FXML lateinit var Pn: TextField
    @FXML lateinit var ro1: TextField
    @FXML lateinit var ro2: TextField
    @FXML lateinit var mu1: TextField
    @FXML lateinit var mu2: TextField
    @FXML lateinit var Gp: TextField
    @FXML lateinit var t: TextField
    @FXML lateinit var n: TextField
    @FXML lateinit var info: TextArea
    @FXML lateinit var tabPane: TabPane
    @FXML fun hydro() {
        val h = Hydro(ro1.text.toDouble(), ro2.text.toDouble(), mu1.text.toDouble(), mu2.text.toDouble(), Gp.text.toDouble(), n.text.toDouble()/100, t.text.toDouble())
        val p = h.calculate(D.text.toDouble()/1000, H1.text.toDouble(), H2.text.toDouble(), L.text.toDouble(), Q.text.toDouble()/24/3600, Pk.text.toDouble()*Math.pow(10.0, 6.0))
        Pn.text = "${(p/Math.pow(10.0, 6.0)*1000).toInt()/1000.0}"
        info.text = h.info
        tabPane.selectionModel.select(2)
    }
}