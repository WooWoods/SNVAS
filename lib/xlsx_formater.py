"""
    xlsx formatter module
    ~~~~~~~~~~~~~~~~~~~~~

    Implements cell format for xlxswriter worksheet.
"""


class Formater:
    def __init__(self, workbook):
        self.header = workbook.add_format(
                {
                    'font_name': 'Times New Roman',
                    'font_size': 12,
                    'align':     'center',
                    'valign':    'vcenter',
                    'border':    True,
                    'bg_color':  'yellow',
                    }
                )

        self.normal = workbook.add_format(
                {
                    'font_name': 'Times New Roman',
                    'font_size': 12,
                    'align':     'center',
                    'valign':    'vcenter',
                    'border':    True,
                    }
                )

        self.remarkable = workbook.add_format(
                {
                    'font_name': 'Times New Roman',
                    'font_size': 12,
                    'align':     'center',
                    'valign':    'vcenter',
                    'border':    True,
                    'bg_color':  '#fac090',
                    }
                )

