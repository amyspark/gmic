/** -*- mode: c++ ; c-basic-offset: 2 -*-
 *
 *  @file LanguageSettings.h
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
#ifndef GMIC_QT_LANGUAGESETTINGS_H
#define GMIC_QT_LANGUAGESETTINGS_H

#include <QMap>
#include <QString>

namespace GmicQt
{

class LanguageSettings {
public:
  LanguageSettings() = delete;
  static const QMap<QString, QString> & availableLanguages();
  static QString configuredTranslator();
  static QString systemDefaultAndAvailableLanguageCode();
  static void installTranslators();
  static bool filterTranslationAvailable(const QString & lang);

private:
  static void installTranslator(const QString & qmPath);
  static void installQtTranslator(const QString & lang);
};

} // namespace GmicQt

#endif // GMIC_QT_LANGUAGESETTINGS_H
