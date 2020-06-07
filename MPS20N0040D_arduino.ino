/*
Pressure Measurements with the
MPS20N0040D Breakout Board
with the HX710B/HX711 ADC
5V Supply Voltage
 */
#include <Q2HX711.h>

const byte MPS_OUT_pin = 2; // OUT data pin
const byte MPS_SCK_pin = 3; // clock data pin
int avg_size = 10; // #pts to average over

Q2HX711 MPS20N0040D(MPS_OUT_pin, MPS_SCK_pin); // start comm with the HX710B

void setup() {
  Serial.begin(9600); // start the serial port
}

void loop() {
  float avg_val = 0.0; // variable for averaging
  for (int ii=0;ii<avg_size;ii++){
    avg_val += MPS20N0040D.read(); // add multiple ADC readings
    delay(50); // delay between readings
  }
  avg_val /= avg_size;
  Serial.println(avg_val,0); // print out the average
}
