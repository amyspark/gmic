/** -*- mode: c++ ; c-basic-offset: 2 -*-
 *
 *  @file CustomDoubleSpinBox.h
 *
 *  Copyright 2017 Sebastien Fourey
 *
 *  This file is part of G'MIC-Qt, a generic plug-in for raster graphics
 *  editors, offering hundreds of filters thanks to the underlying G'MIC
 *  image processing framework.
 *
 *  gmic_qt is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  gmic_qt is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with gmic_qt.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef GMIC_QT_CUSTOMDOUBLESPINBOX_H
#define GMIC_QT_CUSTOMDOUBLESPINBOX_H

#include <QDoubleSpinBox>
#include <QSize>
class QShowEvent;
class QResizeEvent;

class CustomDoubleSpinBox : public QDoubleSpinBox {
  Q_OBJECT
public:
  CustomDoubleSpinBox(QWidget * parent, float min, float max);
  ~CustomDoubleSpinBox() override;
  QString textFromValue(double value) const override;
  inline bool unfinishedKeyboardEditing() const;

protected:
  QSize sizeHint() const override;
  QSize minimumSizeHint() const override;
  void keyPressEvent(QKeyEvent *) override;
signals:
  void keyboardNumericalInputOngoing();

private:
  QSize _sizeHint;
  QSize _minimumSizeHint;
  bool _unfinishedKeyboardEditing = false;
  static const int MAX_DIGITS = 5;
  static int integerPartDigitCount(float value);
};

bool CustomDoubleSpinBox::unfinishedKeyboardEditing() const
{
  return _unfinishedKeyboardEditing;
}

#endif // GMIC_QT_CUSTOMDOUBLESPINBOX_H
