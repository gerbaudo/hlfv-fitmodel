#include <string>
#include <time.h>
#include <stdio.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <iomanip>

using namespace std;

void printTimer(const int eventNumber, 
		const int nrEvents,
		const int modulus)
{
  static time_t begin = time(NULL);
  static time_t now;
  static time_t previous;

  if (eventNumber == 1)
  {
    begin = time(NULL);
  }

  if (eventNumber == 1 || 
      eventNumber % modulus == 0 ||
      eventNumber == nrEvents) {

    previous = now;
    now = time(NULL);
    
    int eventRate = int(begin==now?0.:double(eventNumber)/double(now - begin));
    time_t seconds = (int)eventRate==0?0:(nrEvents-eventNumber)/eventRate;
    time_t minutes = 0;
    time_t hours = 0;
    while (seconds >= 60)
    {
      minutes++;
      seconds -= 60;
    }
    while (minutes >= 60)
    {
      hours++;
      minutes -= 60;
    }



    int instantRate = int(previous==now?0:modulus/(double)(now - previous));
    time_t instantSeconds = (int)instantRate==0?0:(nrEvents-eventNumber)/instantRate;
    time_t instantMinutes = 0;
    time_t instantHours = 0;
    while (instantSeconds >= 60)
    {
      instantMinutes++;
      instantSeconds -= 60;
    }
    while (instantMinutes >= 60)
    {
      instantHours++;
      instantMinutes -= 60;
    }

    cout <<"Event: "<< eventNumber <<"/"<< nrEvents;

    if (eventNumber == 1)
    {
      cout << endl;
    }
    else
    {
      cout << ", Avg/Inst: " << eventRate << "/" << instantRate << " ev/s" 
	   << ", Remaining: " << hours << "h " << minutes << "m " << seconds << "s"
	   << " (" << instantHours << "h " << instantMinutes << "m " << instantSeconds << "s)" << endl;
    }
  }
}
