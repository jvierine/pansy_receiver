
sudo cp sdrreceiver.service /etc/systemd/system/sdrreceiver.service
sudo systemctl daemon-reload
sudo systemctl start sdrreceiver.service
sudo systemctl enable sdrreceiver.service
sudo systemctl status sdrreceiver.service
