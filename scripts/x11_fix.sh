#!/bin/bash
#

# set spaces to pad out subsequent lines printed by "print_" commands
sp11="           "

__log_init__() {
    if [[ -t 1 ]]; then
        # colors for logging in interactive mode
        [[ $COLOR_BOLD ]]   || COLOR_BOLD="\033[1m"
        [[ $COLOR_RED ]]    || COLOR_RED="\033[0;31m"
        [[ $COLOR_GREEN ]]  || COLOR_GREEN="\033[0;34m"
        [[ $COLOR_YELLOW ]] || COLOR_YELLOW="\033[0;33m"
        [[ $COLOR_BLUE ]]   || COLOR_BLUE="\033[0;32m"
        [[ $COLOR_OFF ]]    || COLOR_OFF="\033[0m"
    else
        # no colors to be used if non-interactive
        COLOR_RED= COLOR_GREEN= COLOR_YELLOW= COLOR_BLUE= COLOR_OFF=
    fi
    readonly COLOR_RED COLOR_GREEN COLOR_YELLOW COLOR_BLUE COLOR_OFF

    #
    # map log level strings (FATAL, ERROR, etc.) to numeric values
    #
    # Note the '-g' option passed to declare - it is essential
    #
    unset _log_levels _loggers_level_map
    declare -gA _log_levels _loggers_level_map
    _log_levels=([FATAL]=0 [ERROR]=1 [WARN]=2 [INFO]=3 [DEBUG]=4 [VERBOSE]=5)

}

print_error() {
    {
        printf "${COLOR_RED}ERROR ::   "
        printf '%s\n' "$1"
        [ "$#" -gt "1" ] && { shift; printf "${sp11}%s\n" "$@"; }
        printf "$COLOR_OFF\n"
    } >&2
    exit 1
}

print_warn() {
    printf "${COLOR_YELLOW}WARNING :: "
    printf '%s\n' "$1"
    [ "$#" -gt "1" ] && { shift; printf "${sp11}%s\n" "$@"; }
    printf "$COLOR_OFF"
}

print_info() {
    printf "${COLOR_BLUE}INFO ::    "
    printf '%s\n' "$1"
    [ "$#" -gt "1" ] && { shift; printf "${sp11}%s\n" "$@"; }
    printf "$COLOR_OFF"
}

# print only if output is going to terminal
print_tty() {
    if [[ -t 1 ]]; then
        printf "${sp11}%s\n" "$@"
    fi
}


# get username
username=`who am i | awk '{print $1}'`

# if running as root, add X11 forwarding of this user to root's .Xauthority
if [ "$(id -u)" -eq "0" ]
then
   for (( lineno = 1; lineno <= `xauth -f /home/$username/.Xauthority list | wc -l`; lineno++ ))
   do
      line=`xauth -f /home/$username/.Xauthority list | awk -v ln=$lineno 'NR==ln'`
      xauth add $line
   done
fi

# check X11 forwarding is configured correctly - if not, configure it
if ! xset q &>/dev/null; then
   print_warn "No X server at \$DISPLAY [$DISPLAY]" "Configuring X11 forwarding..."
   xauthstring=`xauth list $DISPLAY`
   xauth add $xauthstring
   export DISPLAY=localhost:10.0
else
   print_info "X11 forwarding configured correctly"
fi

