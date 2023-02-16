/** -*- mode: c++ ; c-basic-offset: 2 -*-
 *
 *  @file IntParameter.cpp
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
#include "FilterParameters/IntParameter.h"
#include <QGridLayout>
#include <QLabel>
#include <QPalette>
#include <QSlider>
#include <QTimerEvent>
#include <QWidget>
#include "FilterParameters/CustomSpinBox.h"
#include "FilterTextTranslator.h"
#include "Globals.h"
#include "HtmlTranslator.h"
#include "Logger.h"
#include "Settings.h"

namespace GmicQt
{

IntParameter::IntParameter(QObject * parent) : AbstractParameter(parent), _min(0), _max(0), _default(0), _value(0), _label(nullptr), _slider(nullptr), _spinBox(nullptr)
{
  _timerId = 0;
  _connected = false;
}

IntParameter::~IntParameter()
{
  delete _spinBox;
  delete _slider;
  delete _label;
}

int IntParameter::size() const
{
  return 1;
}

bool IntParameter::addTo(QWidget * widget, int row)
{
  _grid = dynamic_cast<QGridLayout *>(widget->layout());
  Q_ASSERT_X(_grid, __PRETTY_FUNCTION__, "No grid layout in widget");
  _row = row;
  delete _spinBox;
  delete _slider;
  delete _label;
  _slider = new QSlider(Qt::Horizontal, widget);
  _slider->setMinimumWidth(SLIDER_MIN_WIDTH);
  _slider->setRange(_min, _max);
  _slider->setValue(_value);

  const int delta = 1 + _max - _min;
  if (delta < 20) {
    _slider->setPageStep(1);
  } else {
    const int fact = delta < 100 ? 10 : delta < 1000 ? 100 : delta < 10000 ? 1000 : 10000;
    _slider->setPageStep(fact * (delta / fact) / 10);
  }

  _spinBox = new CustomSpinBox(widget, _min, _max);
  _spinBox->setValue(_value);
#ifndef _GMIC_QT_DISABLE_THEMING_
  if (Settings::darkThemeEnabled()) {
    QPalette p = _slider->palette();
    p.setColor(QPalette::Button, QColor(100, 100, 100));
    p.setColor(QPalette::Highlight, QColor(130, 130, 130));
    _slider->setPalette(p);
  }
#endif
  _grid->addWidget(_label = new QLabel(_name, widget), row, 0, 1, 1);
  setTextSelectable(_label);
  _grid->addWidget(_slider, row, 1, 1, 1);
  _grid->addWidget(_spinBox, row, 2, 1, 1);
  connectSliderSpinBox();
  connect(_spinBox, &CustomSpinBox::editingFinished, [this]() { notifyIfRelevant(); });

  return true;
}

QString IntParameter::value() const
{
  return _spinBox->text();
}

QString IntParameter::defaultValue() const
{
  return QString::number(_default);
}

void IntParameter::setValue(const QString & value)
{
  bool ok = true;
  const int k = value.toInt(&ok);
  if (!ok) {
    Logger::warning(QString("IntParameter::setValue(\"%1\"): bad value").arg(value));
    return;
  }
  _value = k;
  if (_spinBox) {
    disconnectSliderSpinBox();
    _spinBox->setValue(_value);
    _slider->setValue(_value);
    connectSliderSpinBox();
  }
}

void IntParameter::reset()
{
  disconnectSliderSpinBox();
  _slider->setValue(_default);
  _spinBox->setValue(_default);
  _value = _default;
  connectSliderSpinBox();
}

bool IntParameter::initFromText(const QString & filterName, const char * text, int & textLength)
{
  QList<QString> list = parseText("int", text, textLength);
  if (list.isEmpty()) {
    return false;
  }
  _name = HtmlTranslator::html2txt(FilterTextTranslator::translate(list[0], filterName));

  QList<QString> values = list[1].split(QChar(','));
  if (values.size() != 3) {
    return false;
  }
  bool ok1, ok2, ok3;
  _default = values[0].toInt(&ok1);
  _min = values[1].toInt(&ok2);
  _max = values[2].toInt(&ok3);
  _value = _default;
  return ok1 && ok2 && ok3;
}

void IntParameter::timerEvent(QTimerEvent * e)
{
  killTimer(e->timerId());
  _timerId = 0;
  if (!_spinBox->unfinishedKeyboardEditing()) {
    notifyIfRelevant();
  }
}

void IntParameter::onSliderMoved(int value)
{
  if (value != _value) {
    _spinBox->setValue(_value = value);
  }
}

void IntParameter::onSliderValueChanged(int value)
{
  if (value != _value) {
    _spinBox->setValue(_value = value);
  }
}

void IntParameter::onSpinBoxChanged(int i)
{
  _value = i;
  _slider->setValue(i);
  if (_timerId) {
    killTimer(_timerId);
  }
  if (_spinBox->unfinishedKeyboardEditing()) {
    _timerId = 0;
  } else {
    _timerId = startTimer(UPDATE_DELAY);
  }
}

void IntParameter::connectSliderSpinBox()
{
  if (_connected) {
    return;
  }
  connect(_slider, &QSlider::sliderMoved, this, &IntParameter::onSliderMoved);
  connect(_slider, &QSlider::valueChanged, this, &IntParameter::onSliderValueChanged);
  connect(_spinBox, QOverload<int>::of(&CustomSpinBox::valueChanged), this, &IntParameter::onSpinBoxChanged);
  _connected = true;
}

void IntParameter::disconnectSliderSpinBox()
{
  if (!_connected) {
    return;
  }
  _slider->disconnect(this);
  _spinBox->disconnect(this);
  _connected = false;
}

} // namespace GmicQt
