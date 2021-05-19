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

#ifdef Q_OS_ANDROID
#include <android/log.h>
#include <array>
#include <iostream>
#include <thread>
#include <unistd.h>
#endif

#include "DialogSettings.h"
#include "HeadlessProcessor.h"
#include "Host/host.h"
#include "LanguageSettings.h"
#include "Logger.h"
#include "MainWindow.h"
#include "Widgets/ProgressInfoWindow.h"
#include "gmic_qt.h"
#include "gmicqttoolplugin.h"

#include "kpluginfactory.h"

K_PLUGIN_FACTORY_WITH_JSON(KritaGmicPluginFactory, "gmicqttoolplugin.json", registerPlugin<KritaGmicPlugin>();)

KritaGmicPlugin::KritaGmicPlugin(QObject *parent, const QVariantList &)
    : QObject(parent) {}

int KritaGmicPlugin::launch(std::shared_ptr<KisImageInterface> i,
                            bool headless) {
#ifdef Q_OS_ANDROID
  /* Since on Android stdout and stderr redirect to null, un-redirect them */
  /* based on https://stackoverflow.com/a/gmic-qt/31777050 */

  std::array<int, 2> oldFd;
  std::array<int, 2> newStdout, newStderr;

  auto redir_worker = [](std::array<int, 2> &fd, android_LogPriority lvl) {
    ssize_t rdsz;
    std::array<char, 1024> buf{};
    while ((rdsz = read(fd[0], buf.data(), buf.size() - 1)) > 0) {
      if (buf[rdsz - 1] == '\n')
        --rdsz;
      buf[rdsz] = 0; /* add null-terminator */
      __android_log_write(lvl,
                          qPrintable(GmicQt::HostApplicationName), buf.data());
    }
  };

  /* make stdout line-buffered and stderr unbuffered */
  setvbuf(stdout, 0, _IOLBF, 0);
  setvbuf(stderr, 0, _IOLBF, 0);

  /* create the pipe and redirect stdout and stderr */
  dup2(1, oldFd[0]);
  dup2(2, oldFd[1]);
  pipe(newStdout.data());
  pipe(newStderr.data());
  dup2(newStdout[1], 1);
  dup2(newStderr[1], 2);

  /* spawn the logging thread */
  auto newStdoutRedir = std::thread(redir_worker, std::ref(newStdout), ANDROID_LOG_DEBUG);
  auto newStderrRedir = std::thread(redir_worker, std::ref(newStderr), ANDROID_LOG_WARN);
  newStdoutRedir.detach();
  newStderrRedir.detach();
#endif


  // disableInputMode(GmicQt::NoInput);
  // disableInputMode(GmicQt::Active);
  // disableInputMode(GmicQt::All);
  // disableInputMode(GmicQt::ActiveAndBelow);
  // disableInputMode(GmicQt::ActiveAndAbove);
  // disableInputMode(GmicQt::AllVisible);
  // disableInputMode(GmicQt::AllInvisible);

  // disableOutputMode(GmicQt::InPlace);
  disableOutputMode(GmicQt::NewImage);
  disableOutputMode(GmicQt::NewLayers);
  disableOutputMode(GmicQt::NewActiveLayers);

  int r = 0;
  iface = i;
  if (headless) {
    DialogSettings::loadSettings(GmicQt::GuiApplication);
    Logger::setMode(DialogSettings::outputMessageMode());
    LanguageSettings::installTranslators();

    HeadlessProcessor processor;
    QPointer<ProgressInfoWindow> progressWindow = new ProgressInfoWindow(&processor);
    if (processor.command().isEmpty()) {
      return r;
    }

    processor.startProcessing();

    QEventLoop loop;
    connect(progressWindow, SIGNAL(destroyed()), &loop, SLOT(quit()));
    r = loop.exec();

  } else {
    DialogSettings::loadSettings(GmicQt::GuiApplication);
    LanguageSettings::installTranslators();

    QPointer<MainWindow> mainWindow = new MainWindow(0);
    // We want a non modal dialog here.
    mainWindow->setWindowFlags(Qt::Tool | Qt::Dialog);
    mainWindow->setWindowModality(Qt::ApplicationModal);
    mainWindow->setAttribute(Qt::WA_DeleteOnClose); // execution loop never finishes otherwise?

    if (QSettings().value("Config/MainWindowMaximized", false).toBool()) {
      mainWindow->showMaximized();
    }
    else {
      mainWindow->show();
    }

    // Wait than main widget is closed.
    QEventLoop loop;
    connect(mainWindow, SIGNAL(destroyed()), &loop, SLOT(quit()));
    r = loop.exec();
  }

  sharedMemorySegments.clear();
  iface.reset();

#ifdef Q_OS_ANDROID
  /* un-redirect stdout and stderr */
  dup2(oldFd[0], 1);
  dup2(oldFd[1], 2);
#endif

  return r;
}

#include "gmicqttoolplugin.moc"
