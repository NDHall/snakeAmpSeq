#configuration
CONFIG +=  shared qpa no_mocdepend release qt_no_framework
host_build {
    QT_ARCH = x86_64
    QT_TARGET_ARCH = x86_64
} else {
    QT_ARCH = x86_64
    QMAKE_DEFAULT_LIBDIRS = /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/sysroot/lib /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/lib /tmp/build/80754af9/qt_1543467807271/_build_env/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0 /tmp/build/80754af9/qt_1543467807271/_build_env/lib/gcc
    QMAKE_DEFAULT_INCDIRS = /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/include/c++/7.3.0 /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/include/c++/7.3.0/x86_64-conda_cos6-linux-gnu /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/include/c++/7.3.0/backward /tmp/build/80754af9/qt_1543467807271/_build_env/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0/include /tmp/build/80754af9/qt_1543467807271/_build_env/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0/include-fixed /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/include /tmp/build/80754af9/qt_1543467807271/_build_env/x86_64-conda_cos6-linux-gnu/sysroot/usr/include
}
QT_CONFIG +=  minimal-config small-config medium-config large-config full-config fontconfig evdev xlib xrender xcb-plugin xcb-qt xcb-glx xcb-xlib xkbcommon-qt accessibility-atspi-bridge c++11 c++14 c++1z accessibility egl egl_x11 eglfs opengl shared qpa reduce_exports reduce_relocations clock-gettime clock-monotonic posix_fallocate mremap getaddrinfo ipv6ifname getifaddrs inotify eventfd threadsafe-cloexec system-jpeg system-png png system-freetype harfbuzz system-zlib cups iconv glib dbus dbus-linked openssl-linked xcb xinput2 rpath alsa gstreamer-1.0 icu concurrent audio-backend release

#versioning
QT_VERSION = 5.6.3
QT_MAJOR_VERSION = 5
QT_MINOR_VERSION = 6
QT_PATCH_VERSION = 3

#namespaces
QT_LIBINFIX = 
QT_NAMESPACE = 

QT_EDITION = OpenSource

QT_COMPILER_STDCXX = 201402
QT_GCC_MAJOR_VERSION = 7
QT_GCC_MINOR_VERSION = 3
QT_GCC_PATCH_VERSION = 0
