from PyQt5.QtWidgets import (QApplication, QDialog, QComboBox, QLabel, QPushButton,
                             QVBoxLayout, QLineEdit, QMessageBox)
from pydantic import BaseModel, Field, ValidationError, field_validator
from enum import Enum, IntEnum
from typing import List
import sys


# Enum for Color
class ColorEnum(str, Enum):
    Choice1 = 'RED'
    Choice2 = 'GREEN'
    Choice3 = 'BLUE'
    Choice4 = 'CYAN'
    Choice5 = 'Purple'


colormapping = {
    'RED': int,
    'GREEN': float,
    'BLUE': str,
    'CYAN': List[int],
    'Purple': tuple
}


# Pydantic model for validating different input types based on color
class Config(BaseModel):
    color: ColorEnum = Field(title='color', description="Choose a color")
    value: str = Field(title='value', description="Input value based on color")

    # Validator for checking the input type based on the selected color
    @field_validator('value')
    def validate_value_based_on_color(cls, v, values):
        color = values.data['color']
        try:
            if colormapping[color] is List[int]:
                v = [int(i.strip()) for i in v.split(',')]
            elif colormapping[color] is List[float]:
                v = [float(i.strip()) for i in v.split(',')]
            else:
                v = colormapping[color](v)
        except:
            raise ValueError(f'{color} must be {str(colormapping[color])}')
        return v


# PyQt5 dialog with QComboBox and QLineEdit to handle different data types
class ComboBoxWithDynamicValidationApp(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle('Dynamic QComboBox Input Validation')

        # Create a QComboBox
        self.combo = QComboBox()
        # Add Enum values as items
        self.combo.addItems([color.value for color in ColorEnum])
        self.combo.currentIndexChanged.connect(self.on_combo_change)

        # Create QLineEdit for user input
        self.input_field = QLineEdit()
        self.input_field.setPlaceholderText('Enter an integer for Red')

        # Create a QLabel to show instructions
        self.label = QLabel("Select a color and enter corresponding value:")

        # Submit button to validate the selection
        self.submit_button = QPushButton("Submit")
        self.submit_button.clicked.connect(self.validate_selection)

        # Layout setup
        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.combo)
        layout.addWidget(self.input_field)
        layout.addWidget(self.submit_button)

        self.setLayout(layout)

    def on_combo_change(self):
        """Update the input placeholder text based on the selected color."""
        selected_color = self.combo.currentText()

        self.input_field.setPlaceholderText(
            f'Enter {str(colormapping[selected_color])} for {selected_color}')

    def validate_selection(self):
        selected_color = self.combo.currentText()  # Get the current text from QComboBox
        input_value = self.input_field.text()  # Get the input from QLineEdit

        try:
            # Validate the selected item and input value using Pydantic
            config = Config(color=ColorEnum(selected_color), value=input_value)
            QMessageBox.information(
                self, "Success", f"Validated {selected_color}: {input_value}")
        except ValidationError as e:
            QMessageBox.critical(self, "Validation Error", f"Error: {str(e)}")


# Main function to run the PyQt5 app
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = ComboBoxWithDynamicValidationApp()
    window.show()
    sys.exit(app.exec_())
