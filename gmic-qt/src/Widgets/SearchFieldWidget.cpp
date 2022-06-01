/** -*- mode: c++ ; c-basic-offset: 2 -*-
 *
 *  @file SearchFieldWidget.cpp
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
#include "Widgets/SearchFieldWidget.h"
#include <QAction>
#include <QFrame>
#include <QHBoxLayout>
#include <QIcon>
#include <QLineEdit>
#include <QPalette>
#include <QRegExp>
#include <QRegExpValidator>
#include <QToolButton>
#include "Common.h"
#include "IconLoader.h"
#include "Settings.h"
#include "ui_SearchFieldWidget.h"

namespace GmicQt
{

SearchFieldWidget::SearchFieldWidget(QWidget * parent) : QWidget(parent), ui(new Ui::SearchFieldWidget)
{
  ui->setupUi(this);
  _clearIcon = LOAD_ICON("edit-clear");
  _findIcon = LOAD_ICON("edit-find");
  _empty = true;

#if QT_VERSION_GTE(5, 2, 0)
  auto hbox = dynamic_cast<QHBoxLayout *>(layout());
  if (hbox) {
    hbox->setMargin(0);
    hbox->setSpacing(0);
    hbox->addWidget(_lineEdit = new QLineEdit(this));
    _action = _lineEdit->addAction(LOAD_ICON("edit-find"), QLineEdit::TrailingPosition);
    connect(_action, SIGNAL(triggered(bool)), _lineEdit, SLOT(clear()));
  }
#else
  QFrame * frame = new QFrame(this);
  layout()->addWidget(frame);
  frame->setFrameShape(QFrame::StyledPanel);
  frame->setStyleSheet(QString("background:%1").arg(frame->palette().base().color().name()));
  QHBoxLayout * hbox = new QHBoxLayout;
  frame->setLayout(hbox);
  hbox->setMargin(2);

  _lineEdit = new QLineEdit(frame);
  hbox->addWidget(_lineEdit);
  _button = new QToolButton(frame);
  hbox->addWidget(_button);

  _button->setStyleSheet("border:none");
  _lineEdit->setFrame(false);
  _button->setIcon(_findIcon);
  connect(_button, SIGNAL(clicked(bool)), _lineEdit, SLOT(clear()));
#endif
  connect(_lineEdit, SIGNAL(textChanged(QString)), this, SIGNAL(textChanged(QString)));
  connect(_lineEdit, SIGNAL(textChanged(QString)), this, SLOT(onTextChanged(QString)));
  _lineEdit->setPlaceholderText(tr("Search"));
  _lineEdit->setToolTip(tr("Search in filters list (%1)").arg(QKeySequence(QKeySequence::Find).toString()));
  setFocusProxy(_lineEdit);
#if QT_VERSION_GTE(5, 12, 0)
  if (Settings::darkThemeEnabled()) {
    QPalette palette = _lineEdit->palette();
    palette.setColor(QPalette::PlaceholderText, Qt::gray);
    _lineEdit->setPalette(palette);
  }
#endif
  QRegExpValidator * validator = new QRegExpValidator(QRegExp("[^/].*"), this);
  _lineEdit->setValidator(validator);
}

SearchFieldWidget::~SearchFieldWidget()
{
  delete ui;
}

QString SearchFieldWidget::text() const
{
  return _lineEdit->text();
}

void SearchFieldWidget::onTextChanged(const QString & str)
{
#if QT_VERSION_GTE(5, 2, 0)
  if (str.isEmpty()) {
    _empty = true;
    _action->setIcon(_findIcon);
  } else {
    if (_empty) {
      _action->setIcon(_clearIcon);
    }
    _empty = false;
  }
#else
  if (str.isEmpty()) {
    _empty = true;
    _button->setIcon(_findIcon);
  } else {
    if (_empty) {
      _button->setIcon(_clearIcon);
      _empty = false;
    }
  }
#endif
}

void SearchFieldWidget::clear()
{
  _lineEdit->clear();
}

} // namespace GmicQt
