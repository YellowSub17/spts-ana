import sys
from PyQt5.QtWidgets import QApplication, QLabel, QWidget, QVBoxLayout

def run_test():
    # 1. Create the application object
    app = QApplication(sys.argv)

    # 2. Create a basic window
    window = QWidget()
    window.setWindowTitle('X11 Test')
    window.setGeometry(100, 100, 280, 80)

    # 3. Add a "Hello World" label
    layout = QVBoxLayout()
    label = QLabel('<h1>Hello from Qt!</h1>')
    layout.addWidget(label)
    window.setLayout(layout)

    # 4. Show the window and exit
    window.show()
    print("Window should be open. Press Ctrl+C in terminal to exit.")
    sys.exit(app.exec_())

if __name__ == "__main__":
    run_test()
