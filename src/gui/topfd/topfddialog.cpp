//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <map>
#include <string>
#include <vector>
#include <sstream>

#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QDesktopServices>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "base/file_util.hpp"
#include "base/base_data.hpp"
#include "feature/deconv_para.hpp"
#include "feature/deconv_process.hpp"
#include "feature/feature_detect.hpp"

#include "topfddialog.h"
#include "ui_topfddialog.h"
#include "threadtopfd.h"

TopFDDialog::TopFDDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopFDDialog) {
      initArguments();
      ui->setupUi(this);
      lastDir_ = ".";
      percentage_ = 0;
      ui->maxChargeEdit->setValidator(new QIntValidator(1, 100, this));
      ui->maxMassEdit->setValidator(new QIntValidator(1, 1000000, this));
      QRegExp rx1("^0\\.[0]\\d{0,2}[1-9]|0.1$");
      QRegExpValidator *validator1 = new QRegExpValidator(rx1, this);
      ui->mzErrorEdit->setValidator(validator1);
      QRegExp rx2("^\\d{1,6}\\.\\d{0,2}$");
      QRegExpValidator *validator2 = new QRegExpValidator(rx2, this);
      ui->ms1snRatioEdit->setValidator(validator2);
      ui->ms2snRatioEdit->setValidator(validator2);
      QRegExp rx3("^\\d{1,4}\\.\\d{0,2}|10000$");
      QRegExpValidator *validator3 = new QRegExpValidator(rx3, this);
      ui->windowSizeEdit->setValidator(validator3);
      QFont font;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
      font.setFamily(QStringLiteral("Courier New"));
#else
      font.setFamily(QStringLiteral("Monospace"));
#endif
      ui->outputTextBrowser->setFont(font);
      thread_ = new ThreadTopFD(this);
      showInfo = "";
      TopFDDialog::on_defaultButton_clicked();
    }

TopFDDialog::~TopFDDialog() {
  thread_->terminate();
  thread_->wait();
  delete ui;
}

void TopFDDialog::closeEvent(QCloseEvent *event) {
  if (thread_->isRunning()) {
    if (!continueToClose()) {
      event->ignore();
      return;
    }
  }
  thread_->terminate();
  thread_->wait();
  event->accept();
  return;
}

void TopFDDialog::initArguments() {
  arguments_["executiveDir"] = "";
  arguments_["resourceDir"] = "";
  arguments_["spectrumFileName"] = "";
  arguments_["refinePrecMass"]="true";
  arguments_["missingLevelOne"] = "false";
  arguments_["maxCharge"] = "30";
  arguments_["maxMass"] = "100000";
  arguments_["mzError"] = "0.02";
  arguments_["msTwoSnRatio"] = "1.0";
  arguments_["msOneSnRatio"] = "3.0";
  arguments_["keepUnusedPeaks"] = "false";
  arguments_["outMultipleMass"] = "false";
  arguments_["precWindow"] = "3.0";
  arguments_["doFinalFiltering"] = "true";
  arguments_["outputMatchEnv"] = "false";
}

void TopFDDialog::on_clearButton_clicked() {
  //ui->maxChargeEdit->clear();
  //ui->maxMassEdit->clear();
  //ui->mzErrorEdit->clear();
  //ui->ms1snRatioEdit->clear();
  //ui->ms2snRatioEdit->clear();
  //ui->windowSizeEdit->clear();
  //ui->missLevelOneCheckBox->setChecked(false);
  ui->listWidget->clear();
  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  lastDir_ = "/";
}

void TopFDDialog::on_defaultButton_clicked() {
  ui->maxChargeEdit->setText("30");
  ui->maxMassEdit->setText("100000");
  ui->mzErrorEdit->setText("0.02");
  ui->ms1snRatioEdit->setText("3.0");
  ui->ms2snRatioEdit->setText("1.0");
  ui->windowSizeEdit->setText("3.0");
  ui->missLevelOneCheckBox->setChecked(false);
  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
}

std::vector<std::string> TopFDDialog::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void TopFDDialog::on_addButton_clicked() {
  QStringList spfiles = QFileDialog::getOpenFileNames(
      this,
      "Select spectrum files",
      lastDir_,
      "Spectra files(*.mzXML *.mzML *.mzxml *.mzml)");
  for (int i = 0; i < spfiles.size(); i++) {
    QString spfile = spfiles.at(i);
    updatedir(spfile);
    if (ableToAdd(spfile)) {
      ui->listWidget->addItem(new QListWidgetItem(spfile));
    }
  }
}

void TopFDDialog::updatedir(QString s) {
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
}
bool TopFDDialog::ableToAdd(QString spfile) {
  bool able = true;
  if (spfile != "") {
    if (spfile.toStdString().length() > 200) {
      QMessageBox::warning(this, tr("Warning"),
                           tr("The spectrum file path is too long!"),
                           QMessageBox::Yes);
      able = false;
    } else {
      for (int i = 0; i < ui->listWidget->count(); i++) {
        if (spfile == ui->listWidget->item(i)->text()) {
          able = false;
        }
      }
    }
  } else {
    able = false;
  }
  return able;
}

void TopFDDialog::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
}

void TopFDDialog::on_startButton_clicked() {
  std::stringstream buffer;
  std::streambuf *oldbuf = std::cout.rdbuf(buffer.rdbuf());
  if (checkError()) {
    return;
  }
  lockDialog();
  ui->outputTextBrowser->setText(showInfo);
  std::map<std::string, std::string> argument = this->getArguments();
  std::vector<std::string> spec_file_lst = this->getSpecFileList();
  thread_->setPar(argument, spec_file_lst);
  thread_->start();

  std::string info;
  int processed_len = 0;
  std::string processed_lines = ""; 
  std::string current_line = "";
  unsigned cursor_pos = 0;

  while (true) {
    // Here is the infomation been shown in the infoBox.
    info = buffer.str();
    std::string new_info = info.substr(processed_len);
    processed_len = info.length();
    
    if (new_info.size() > 0) {
      for (unsigned i = 0; i < new_info.size(); i++) {
        // new line
        if (new_info.at(i) == '\n') {
          processed_lines = processed_lines + current_line + '\n';
          current_line = "";
          cursor_pos = 0;
        }
        // CF
        if (new_info.at(i) == '\r') {
          cursor_pos = 0;
        }
        // add a new charactor
        if (new_info.at(i) != '\n' && new_info.at(i) != '\r') {
          if (cursor_pos < current_line.length()) {
            current_line[cursor_pos] = new_info.at(i);
          }
          else {
            current_line = current_line + new_info.at(i);
          }
          cursor_pos++;
        }
      }
      updateMsg(processed_lines + current_line);
    }
    if (thread_->isFinished()) {
      break;
    }
    sleep(100);
  }
  unlockDialog();
  showInfo = "";
  thread_->exit();
  std::cout.rdbuf(oldbuf);
}

void TopFDDialog::on_exitButton_clicked() {
  close();
}

bool TopFDDialog::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopFD is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}

void TopFDDialog::on_outputButton_clicked() {
  std::string sp_file_name = "";
  if (spec_file_lst_.size() > 0) {
    sp_file_name = spec_file_lst_[0];
  }
  fs::path full_path(sp_file_name.c_str());
  QString outPath = full_path.remove_filename().string().c_str();
  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

std::map<std::string, std::string> TopFDDialog::getArguments() {
  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = prot::file_util::getExecutiveDir(path.toStdString());
  arguments_["executiveDir"] = exe_dir;
  arguments_["resourceDir"] = arguments_["executiveDir"] + prot::file_util::getFileSeparator() + prot::file_util::getResourceDirName();
  arguments_["maxCharge"] = ui->maxChargeEdit->text().toStdString();
  arguments_["maxMass"] = ui->maxMassEdit->text().toStdString();
  arguments_["mzError"] = ui->mzErrorEdit->text().toStdString();
  arguments_["msOneSnRatio"] = ui->ms1snRatioEdit->text().toStdString();
  arguments_["msTwoSnRatio"] = ui->ms2snRatioEdit->text().toStdString();
  arguments_["precWindow"] = ui->windowSizeEdit->text().toStdString();
  return arguments_;
}

void TopFDDialog::lockDialog() {
  ui->addButton->setEnabled(false);
  ui->delButton->setEnabled(false);
  ui->maxChargeEdit->setEnabled(false);
  ui->maxMassEdit->setEnabled(false);
  ui->mzErrorEdit->setEnabled(false);
  ui->ms1snRatioEdit->setEnabled(false);
  ui->ms2snRatioEdit->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->missLevelOneCheckBox->setEnabled(false);
  ui->windowSizeEdit->setEnabled(false);
  ui->outputButton->setEnabled(false);
}

void TopFDDialog::unlockDialog() {
  ui->addButton->setEnabled(true);
  ui->delButton->setEnabled(true);
  ui->maxChargeEdit->setEnabled(true);
  ui->maxMassEdit->setEnabled(true);
  ui->mzErrorEdit->setEnabled(true);
  ui->ms1snRatioEdit->setEnabled(true);
  ui->ms2snRatioEdit->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->missLevelOneCheckBox->setEnabled(true);
  ui->windowSizeEdit->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
}

bool TopFDDialog::checkError() {
  if (ui->maxChargeEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Maximum charge is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->maxMassEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Maximum mass is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->mzErrorEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("M/z error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms1snRatioEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS1 S/N ratio is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms2snRatioEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS2 S/N ratio is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->windowSizeEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Precursor window size is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}

void TopFDDialog::updateMsg(std::string msg) {
  showInfo = msg.c_str();
  ui->outputTextBrowser->setText(showInfo);
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
}

void TopFDDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}

