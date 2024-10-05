from PyQt5.QtWidgets import (QApplication, QDialog, QGridLayout, QLabel, QComboBox,
                             QLineEdit, QRadioButton, QHBoxLayout, QPushButton, QVBoxLayout, QMessageBox,QFileDialog)
from pydantic import BaseModel, Field, ValidationError, DirectoryPath, FilePath
from enum import IntEnum
from typing import Optional, get_origin, Union, List, get_args
import sys
from pathlib import Path
import os


class IntEnumWithAlias(IntEnum):
    def __new__(cls, value, alias):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.alias = alias
        return obj


class ItemEnum(IntEnumWithAlias):
    APPLE = 0, 'APPLE'
    BANANA = 1, 'BANANA'
    ORANGE = 2, 'ORANGE'


class Config(BaseModel):
    workspace: DirectoryPath =Field(alias='workspace',description="workspace",default=f'{os.getcwd()}')
    Address: str = Field(alias='Address', default='taiwan',
                         description='Address')
    port: int = Field(alias='door', default=56, description='door_num')
    Option1: Optional[float] = Field(
        alias="Sour", description='Sour', default=None)
    Option2: Optional[int] = Field(
        alias='Salinity', description='Salinity', default=None)
    Option3: Optional[List[int]] = Field(
        alias='Sweetness', description='Sweetness', default=None)
    item: ItemEnum = Field(
        alias='item', default=ItemEnum.APPLE, description='components')
    f: FilePath = Field(alias='FilePath', description='structure_path',default=None)


class PydanticDialog(QDialog):
    def __init__(self, model: type[BaseModel], parent=None):
        super().__init__(parent=parent)
        self.model = model
        self.pydantic_layout = QGridLayout()
        self.pydantic_widgets = dict()
        self.constructFromMetaInfo()

    def constructFromMetaInfo(self):
        fields = self.model.model_fields
        for row, (field_name, field_info) in enumerate(fields.items()):
            self.pydantic_layout.addWidget(
                QLabel(field_info.description), row, 0)
            if get_origin(field_info.annotation) is None:
                if issubclass(field_info.annotation, IntEnum):
                    self.pydantic_widgets[field_name] = QComboBox()
                    self.pydantic_widgets[field_name].addItems(
                        [item.alias for item in list(field_info.annotation)])
                    default_value = field_info.default.alias if field_info.default else ''
                    self.pydantic_widgets[field_name].setCurrentText(
                        default_value)
                    self.pydantic_layout.addWidget(
                        self.pydantic_widgets[field_name], row, 1)
                    
                elif issubclass(field_info.annotation,Path):
                    self.pydantic_widgets[field_name]=QLineEdit()
                    file_button = QPushButton("Browse")
                    self.pydantic_layout.addWidget(file_button,row,2)
                    self.pydantic_layout.addWidget(self.pydantic_widgets[field_name],row,1)
                    pathtype=field_info.metadata[0].path_type
                    file_button.clicked.connect(lambda _,x=field_name,y=pathtype:self.browseFile(x,y))

                elif issubclass(field_info.annotation, str):
                    self.pydantic_widgets[field_name] = QLineEdit()
                    self.pydantic_layout.addWidget(
                        self.pydantic_widgets[field_name], row, 1)
                elif issubclass(field_info.annotation, float):
                    self.pydantic_widgets[field_name] = QLineEdit()
                    self.pydantic_layout.addWidget(
                        self.pydantic_widgets[field_name], row, 1)
                elif issubclass(field_info.annotation, int):
                    self.pydantic_widgets[field_name] = QLineEdit()
                    self.pydantic_layout.addWidget(
                        self.pydantic_widgets[field_name], row, 1)

            elif get_origin(field_info.annotation) is Union:
                # Adding a QRadioButton for optional fields
                if 'Option' in field_name:
                    layout = QHBoxLayout()
                    radio_button = QRadioButton("Enable")
                    edit_line = QLineEdit()
                    edit_line.setEnabled(False)

                    # Connect radio button toggling to enable/disable the edit line
                    radio_button.toggled.connect(
                        lambda checked, edit_line=edit_line: edit_line.setEnabled(checked))

                    layout.addWidget(radio_button)
                    layout.addWidget(edit_line)
                    self.pydantic_layout.addLayout(layout, row, 1)
                    # Store the radio button and line edit as part of the widget for later retrieval
                    self.pydantic_widgets[field_name] = (
                        radio_button, edit_line)
                else:
                    self.pydantic_widgets[field_name] = QLineEdit()
                    self.pydantic_layout.addWidget(
                        self.pydantic_widgets[field_name], row, 1)

    def asLayout(self):
        self.setLayout(self.pydantic_layout)

    def browseFile(self,field_name,pathtype):
        file_dialog=QFileDialog()
        if pathtype=='file':
            file_dialog.setFileMode(QFileDialog.AnyFile)
            file_dialog.setAcceptMode(QFileDialog.AcceptOpen)
            file_path,_=file_dialog.getOpenFileName(self,'Select File')
        elif pathtype=='dir':
            file_dialog.setFileMode(QFileDialog.DirectoryOnly)
            file_dialog.setAcceptMode(QFileDialog.AcceptOpen)
            file_path=file_dialog.getExistingDirectory(self,'Select Directory')
        else:
             print(f"Unknown path type: {pathtype}")
             return 
        if file_path:
            self.pydantic_widgets[field_name].setText(file_path)

    def getData(self):
        data = {}
        for field_name, widget in self.pydantic_widgets.items():
            # for optional fields (with QRadioButton and QLineEdit)
            if isinstance(widget, tuple):
                radio_button, edit_line = widget
                if radio_button.isChecked():
                    value = edit_line.text()
                    # We have to tackle with list
                    field_annotation = self.model.model_fields[field_name].annotation
                    if get_args(field_annotation)[0] is List[int]:
                        data[field_name] = [int(i.strip())
                                            for i in value.split(',') if i.strip()]
                    elif get_args(field_annotation)[0] is List[float]:
                        data[field_name] = [float(i.strip())
                                            for i in value.split(',') if i.strip()]
                    else:
                        data[field_name] = value
            elif isinstance(widget, QComboBox):
                selected_item = widget.currentText()
                data[field_name] = ItemEnum[selected_item].value
            elif isinstance(widget, QLineEdit):
                data[field_name] = widget.text()

        return {self.model.model_fields[field].alias: value for field, value in data.items()}


class ConfigDialog(PydanticDialog):
    def __init__(self, parent=None):
        super().__init__(Config, parent)
        self.setWindowTitle('Pydantic_PyQt5')

        # Add QPushButton
        self.submit_button = QPushButton('Submit')
        self.submit_button.clicked.connect(self.handle_submit)
        # Create a vertical layout to include the form and the button
        main_layout = QVBoxLayout()
        main_layout.addLayout(self.pydantic_layout)
        main_layout.addWidget(self.submit_button)

        self.setLayout(main_layout)

    def handle_submit(self):
        data = self.getData()
        try:
            # Validate data using base_model
            config = Config.model_validate(data)
            print(f'Validated Config Object:{config}')
            QMessageBox.information(
                self, 'Success', 'Data has been successfully validated!')
        except ValidationError as e:
            QMessageBox.critical(self, 'Validation Error', str(e))
            print(f'Validation Error {e}')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    cd = ConfigDialog()
    cd.show()
    sys.exit(app.exec_())
