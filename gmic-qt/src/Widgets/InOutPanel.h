/** -*- mode: c++ ; c-basic-offset: 2 -*-
 *
 *  @file InOutPanel.h
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
#ifndef GMIC_QT_INOUTPANEL_H
#define GMIC_QT_INOUTPANEL_H

#include <QWidget>
#include "Host/GmicQtHost.h"
#include "InputOutputState.h"
#include "GmicQt.h"
class QSettings;
class QPalette;

namespace Ui
{
class InOutPanel;
}

namespace GmicQt
{

class FilterThread;

class InOutPanel : public QWidget {
  Q_OBJECT

public:
  explicit InOutPanel(QWidget * parent = nullptr);
  ~InOutPanel();

public:
  InputMode inputMode() const;
  OutputMode outputMode() const;
  OutputMessageMode outputMessageMode() const;
  void reset();

  void disableNotifications();
  void enableNotifications();
  void setInputMode(InputMode mode);
  void setOutputMode(OutputMode mode);

  InputOutputState state() const;
  void setState(const InputOutputState & state, bool notify);

  void setEnabled(bool);
  void disable();
  void enable();

  static void disableInputMode(InputMode mode);
  static void disableOutputMode(OutputMode mode);

  bool hasActiveControls();

signals:
  void inputModeChanged(InputMode);

public slots:
  void onInputModeSelected(int);
  void onOutputModeSelected(int);
  void onResetButtonClicked();
#ifndef _GMIC_QT_DISABLE_THEMING_
  void setDarkTheme();
#endif

private:
  static void setDefaultInputMode();
  static void setDefaultOutputMode();
  void setTopLabel();
  void updateLayoutIfUniqueRow();
  bool _notifyValueChange;
  Ui::InOutPanel * ui;
  static const int NoSelection = -1;
  static QList<InputMode> _enabledInputModes;
  static QList<OutputMode> _enabledOutputModes;
};

} // namespace GmicQt

#endif // GMIC_QT_INOUTPANEL_H
