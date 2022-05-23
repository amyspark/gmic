.PHONY: all clean

QMS = de.qm es.qm fr.qm it.qm ja.qm nl.qm pl.qm pt.qm ru.qm zh.qm

.PRECIOUS: $(patsubst %.qm,%.ts,$(QMS))

all: $(QMS)

%.qm: %.ts
	@../lrelease.sh $<

%.ts: gmic_qt_%.csv
	@echo Generating $@ file from $<
	@./csv2ts.sh -o $@ $<

clean:
	rm -f $(QMS) *.ts
