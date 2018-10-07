import javafx.application.Application
import javafx.scene.Scene
import javafx.stage.Stage
import javafx.fxml.FXMLLoader
import javafx.scene.Parent



@Suppress("JAVA_CLASS_ON_COMPANION")
class Main : Application() {
    override fun start(primaryStage: Stage) {
        primaryStage.title = "Гидравлический расчёт"
        val root = FXMLLoader.load<Parent>(javaClass.getResource("form.fxml"))
        primaryStage.scene = Scene(root)
        primaryStage.show()
    }
    companion object {
        @JvmStatic
        fun main(args: Array<String>) {
            launch(Main::class.java)
        }
    }
}
