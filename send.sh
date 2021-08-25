PASSWORD="xxxxxxxx"

FILENAME="solve.cpp"
#FILENAME="run.sh"

expect -c "
set timeout 10
spawn scp ${FILENAME} l00164@login.fugaku.r-ccs.riken.jp:${FILENAME}
expect \"Enter passphrase for key '/Users/sk/.ssh/id_ed25519':\"
send \"${PASSWORD}\n\"
expect \"$\"
"
