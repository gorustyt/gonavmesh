package gui

const (
	MAX_MESSAGES   = 1000
	TEXT_POOL_SIZE = 8000
)

type logger struct {
	m_messageCount int
	m_textPoolSize int
	m_messages     [MAX_MESSAGES]string
}

func newLogger() *logger {
	return &logger{}
}

// / Dumps the log to stdout.
func (l *logger) dumpLog(format string) {

}

// / Returns number of log messages.
func (l *logger) getLogCount() int {
	return l.m_messageCount
}

// / Returns log message text.
func (l *logger) getLogText(i int) string {
	return l.m_messages[i]
}

func (l *logger) reset() {
	l.m_messageCount = 0
}
