mkdir -p ~/.config/systemd/user/

cp sdrreceiver.service ~/.config/systemd/user/sdrreceiver.service
systemctl --user daemon-reload
systemctl --user start sdrreceiver.service
systemctl --user enable sdrreceiver.service
systemctl --user status sdrreceiver.service

