/** -*- mode: c++ ; c-basic-offset: 2 -*-
 *
 *  @file Settings.cpp
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
#include "Settings.h"
#include "Host/GmicQtHost.h"
#include "IconLoader.h"
namespace
{
GmicQt::OutputMessageMode filterDeprecatedOutputMessageMode(const GmicQt::OutputMessageMode & mode)
{
  if (mode == GmicQt::OutputMessageMode::VerboseLayerName_DEPRECATED) {
    return GmicQt::DefaultOutputMessageMode;
  }
  return mode;
}
} // namespace

namespace GmicQt
{

bool Settings::_visibleLogos;
bool Settings::_darkThemeEnabled;
QString Settings::_languageCode;
bool Settings::_filterTranslationEnabled = false;
MainWindow::PreviewPosition Settings::_previewPosition;
bool Settings::_nativeColorDialogs;
int Settings::_updatePeriodicity;
int Settings::_previewTimeout = 16;
OutputMessageMode Settings::_outputMessageMode;
bool Settings::_previewZoomAlwaysEnabled = false;
bool Settings::_notifyFailedStartupUpdate = true;

const QColor Settings::CheckBoxBaseColor(83, 83, 83);
const QColor Settings::CheckBoxTextColor(255, 255, 255);
QColor Settings::UnselectedFilterTextColor;

QString Settings::FolderParameterDefaultValue;
QString Settings::FileParameterDefaultPath;

QIcon Settings::AddIcon;
QIcon Settings::RemoveIcon;
QString Settings::GroupSeparator;
QString Settings::DecimalPoint('.');
QString Settings::NegativeSign('-');

void Settings::load(UserInterfaceMode userInterfaceMode)
{
  QSettings settings;
  _visibleLogos = settings.value("LogosAreVisible", true).toBool();
  _darkThemeEnabled = settings.value(DARK_THEME_KEY, GmicQtHost::DarkThemeIsDefault).toBool();
  _languageCode = settings.value(LANGUAGE_CODE_KEY, QString()).toString();

  if (settings.value("Config/PreviewPosition", "Left").toString() == "Left") {
    _previewPosition = MainWindow::PreviewPosition::Left;
  } else {
    _previewPosition = MainWindow::PreviewPosition::Right;
  }
  _filterTranslationEnabled = settings.value(ENABLE_FILTER_TRANSLATION, false).toBool();
  _nativeColorDialogs = settings.value("Config/NativeColorDialogs", false).toBool();
  _updatePeriodicity = settings.value(INTERNET_UPDATE_PERIODICITY_KEY, INTERNET_DEFAULT_PERIODICITY).toInt();
  FolderParameterDefaultValue = settings.value("FolderParameterDefaultValue", QDir::homePath()).toString();
  FileParameterDefaultPath = settings.value("FileParameterDefaultPath", QDir::homePath()).toString();
  _previewTimeout = settings.value("PreviewTimeout", 16).toInt();
  _previewZoomAlwaysEnabled = settings.value("AlwaysEnablePreviewZoom", false).toBool();
  _outputMessageMode = filterDeprecatedOutputMessageMode((GmicQt::OutputMessageMode)settings.value("OutputMessageMode", static_cast<int>(GmicQt::DefaultOutputMessageMode)).toInt());
  _notifyFailedStartupUpdate = settings.value("Config/NotifyIfStartupUpdateFails", true).toBool();
  if (userInterfaceMode != UserInterfaceMode::Silent) {
    AddIcon = LOAD_ICON("list-add");
    RemoveIcon = LOAD_ICON("list-remove");
  }
  QLocale locale;
  GroupSeparator = locale.groupSeparator();
  DecimalPoint = locale.decimalPoint();
  NegativeSign = locale.negativeSign();
}

bool Settings::visibleLogos()
{
  return _visibleLogos;
}

void Settings::setVisibleLogos(bool on)
{
  _visibleLogos = on;
}

bool Settings::darkThemeEnabled()
{
#ifdef _GMIC_QT_DISABLE_THEMING_
  return GmicQtHost::DarkThemeIsDefault;
#else
  return _darkThemeEnabled;
#endif
}

void Settings::setDarkThemeEnabled(bool on)
{
  _darkThemeEnabled = on;
}

QString Settings::languageCode()
{
  return _languageCode;
}

void Settings::setLanguageCode(const QString & code)
{
  _languageCode = code;
}

bool Settings::filterTranslationEnabled()
{
  return _filterTranslationEnabled;
}

void Settings::setFilterTranslationEnabled(bool on)
{
  _filterTranslationEnabled = on;
}

MainWindow::PreviewPosition Settings::previewPosition()
{
  return _previewPosition;
}

void Settings::setPreviewPosition(MainWindow::PreviewPosition position)
{
  _previewPosition = position;
}

bool Settings::nativeColorDialogs()
{
  return _nativeColorDialogs;
}

void Settings::setNativeColorDialogs(bool on)
{
  _nativeColorDialogs = on;
}

int Settings::updatePeriodicity()
{
  return _updatePeriodicity;
}

void Settings::setUpdatePeriodicity(int hours)
{
  _updatePeriodicity = hours;
}

int Settings::previewTimeout()
{
  return _previewTimeout;
}

void Settings::setPreviewTimeout(int seconds)
{
  _previewTimeout = seconds;
}

OutputMessageMode Settings::outputMessageMode()
{
  return _outputMessageMode;
}

void Settings::setOutputMessageMode(OutputMessageMode mode)
{
  _outputMessageMode = mode;
}

bool Settings::previewZoomAlwaysEnabled()
{
  return _previewZoomAlwaysEnabled;
}

void Settings::setPreviewZoomAlwaysEnabled(bool on)
{
  _previewZoomAlwaysEnabled = on;
}

bool Settings::notifyFailedStartupUpdate()
{
  return _notifyFailedStartupUpdate;
}

void Settings::setNotifyFailedStartupUpdate(bool on)
{
  _notifyFailedStartupUpdate = on;
}

void Settings::save(QSettings & settings)
{
  removeObsoleteKeys(settings);
  settings.setValue("LogosAreVisible", _visibleLogos);
  settings.setValue(DARK_THEME_KEY, _darkThemeEnabled);
  settings.setValue(LANGUAGE_CODE_KEY, _languageCode);
  settings.setValue(ENABLE_FILTER_TRANSLATION, _filterTranslationEnabled);
  settings.setValue("Config/PreviewPosition", (_previewPosition == MainWindow::PreviewPosition::Left) ? "Left" : "Right");

  settings.setValue("Config/NativeColorDialogs", _nativeColorDialogs);
  settings.setValue(INTERNET_UPDATE_PERIODICITY_KEY, _updatePeriodicity);
  settings.setValue("FolderParameterDefaultValue", FolderParameterDefaultValue);
  settings.setValue("FileParameterDefaultPath", FileParameterDefaultPath);
  settings.setValue("PreviewTimeout", _previewTimeout);
  settings.setValue("OutputMessageMode", (int)_outputMessageMode);
  settings.setValue("AlwaysEnablePreviewZoom", _previewZoomAlwaysEnabled);
  // Remove obsolete keys (2.0.0 pre-release)
  settings.remove("Config/UseFaveInputMode");
  settings.remove("Config/UseFaveOutputMode");
  settings.remove("Config/UseFaveOutputMessages");
  settings.remove("Config/UseFavePreviewMode");
  settings.setValue("Config/NotifyIfStartupUpdateFails", _notifyFailedStartupUpdate);
}

void Settings::removeObsoleteKeys(QSettings & settings)
{
  settings.remove(QString("LastExecution/host_%1/PreviewMode").arg(GmicQtHost::ApplicationShortname));
  settings.remove(QString("LastExecution/host_%1/GmicEnvironment").arg(GmicQtHost::ApplicationShortname));
  settings.remove(QString("LastExecution/host_%1/QuotedParameters").arg(GmicQtHost::ApplicationShortname));
  settings.remove(QString("LastExecution/host_%1/GmicStatus").arg(GmicQtHost::ApplicationShortname));
}

} // namespace GmicQt
