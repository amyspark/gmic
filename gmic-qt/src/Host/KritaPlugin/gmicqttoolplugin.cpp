/*
 *  This file is part of G'MIC-Qt, a generic plug-in for raster graphics
 *  editors, offering hundreds of filters thanks to the underlying G'MIC
 *  image processing framework.
 *
 *  Copyright (C) 2020 L. E. Segovia <amy@amyspark.me>
 *
 *  Description: Krita painting suite plugin for G'Mic-Qt.
 *
 *  G'MIC-Qt is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  G'MIC-Qt is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <QApplication>
#include <QEventLoop>
#include <QPointer>
#include <QSettings>
#include <QTranslator>
#include <list>

#include "Settings.h"
#include "GmicQt.h"
#include "HeadlessProcessor.h"
#include "Host/GmicQtHost.h"
#include "LanguageSettings.h"
#include "Logger.h"
#include "MainWindow.h"
#include "Widgets/InOutPanel.h"
#include "Widgets/ProgressInfoWindow.h"
#include "gmicqttoolplugin.h"

#include "kpluginfactory.h"

K_PLUGIN_FACTORY_WITH_JSON(KritaGmicPluginFactory,
                           "gmicqttoolplugin.json",
                           registerPlugin<KritaGmicPlugin>();)

KritaGmicPlugin::KritaGmicPlugin(QObject *parent, const QVariantList &)
    : QObject(parent)
{
}

int KritaGmicPlugin::launch(std::shared_ptr<KisImageInterface> i, bool headless)
{
  using namespace GmicQt;

  std::list<GmicQt::InputMode> disabledInputModes;
  disabledInputModes.push_back(GmicQt::InputMode::NoInput);
  // disabledInputModes.push_back(GmicQt::Active);
  // disabledInputModes.push_back(GmicQt::All);
  // disabledInputModes.push_back(GmicQt::ActiveAndBelow);
  // disabledInputModes.push_back(GmicQt::ActiveAndAbove);
  disabledInputModes.push_back(GmicQt::InputMode::AllVisible);
  disabledInputModes.push_back(GmicQt::InputMode::AllInvisible);

  std::list<GmicQt::OutputMode> disabledOutputModes;
  // disabledOutputModes.push_back(GmicQt::OutputMode::InPlace);
  disabledOutputModes.push_back(GmicQt::OutputMode::NewImage);
  disabledOutputModes.push_back(GmicQt::OutputMode::NewLayers);
  disabledOutputModes.push_back(GmicQt::OutputMode::NewActiveLayers);

  int status = 0;
  GmicQtHost::iface = i;
  if (headless) {
    GmicQt::RunParameters parameters = GmicQt::lastAppliedFilterRunParameters(
        GmicQt::ReturnedRunParametersFlag::AfterFilterExecution);
    {
      for (const GmicQt::InputMode & mode : disabledInputModes) {
        GmicQt::InOutPanel::disableInputMode(mode);
      }
      for (const GmicQt::OutputMode & mode : disabledOutputModes) {
        GmicQt::InOutPanel::disableOutputMode(mode);
      }
    }
    Settings::load(GmicQt::UserInterfaceMode::ProgressDialog);
    Logger::setMode(Settings::outputMessageMode());
    LanguageSettings::installTranslators();

    HeadlessProcessor processor(nullptr);
    if (!processor.setPluginParameters(parameters)) {
      Logger::error(processor.error());
      return 1;
    }

    QPointer<ProgressInfoWindow> progressWindow(new ProgressInfoWindow(&processor));
    // We want a non modal dialog here.
    progressWindow->setWindowFlags(Qt::Tool | Qt::Dialog);
    progressWindow->setWindowModality(Qt::ApplicationModal);
    // Make it destroy itself on close (signaling the event loop)
    progressWindow->setAttribute(Qt::WA_DeleteOnClose);

    processor.startProcessing();

    QEventLoop loop;
    connect(progressWindow, SIGNAL(destroyed()), &loop, SLOT(quit()));
    status = loop.exec();
  } else {
    GmicQt::RunParameters parameters = GmicQt::lastAppliedFilterRunParameters(
        GmicQt::ReturnedRunParametersFlag::AfterFilterExecution);
    {
      for (const GmicQt::InputMode & mode : disabledInputModes) {
        GmicQt::InOutPanel::disableInputMode(mode);
      }
      for (const GmicQt::OutputMode & mode : disabledOutputModes) {
        GmicQt::InOutPanel::disableOutputMode(mode);
      }
    }
    Settings::load(GmicQt::UserInterfaceMode::Full);
    Logger::setMode(Settings::outputMessageMode());
    LanguageSettings::installTranslators();

    QPointer<MainWindow> mainWindow(new MainWindow(qApp->activeWindow()));
    mainWindow->setPluginParameters(parameters);
#ifdef Q_OS_MACOS
    mainWindow->setWindowFlags(Qt::Tool | Qt::Dialog);
#else
    mainWindow->setWindowFlags(Qt::Dialog);
#endif
    mainWindow->setWindowModality(Qt::ApplicationModal);
    // Make it destroy itself on close (signaling the event loop)
    mainWindow->setAttribute(Qt::WA_DeleteOnClose);

    if (QSettings().value("Config/MainWindowMaximized", false).toBool()) {
      mainWindow->showMaximized();
    } else {
      mainWindow->show();
    }

    // Wait than main widget is closed.
    QEventLoop loop;
    connect(mainWindow, SIGNAL(destroyed()), &loop, SLOT(quit()));
    status = loop.exec();
  }

  GmicQtHost::sharedMemorySegments.clear();
  GmicQtHost::iface.reset();

  return status;
}

#include "gmicqttoolplugin.moc"
