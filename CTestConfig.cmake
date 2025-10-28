set(CTEST_PROJECT_NAME Pele)
set(CTEST_NIGHTLY_START_TIME 20:00:00 America/Denver)

if(CMAKE_VERSION VERSION_GREATER 3.14)
  set(CTEST_SUBMIT_URL https://my.cdash.org/submit.php?project=Pele)
else()
  set(CTEST_DROP_METHOD "https")
  set(CTEST_DROP_SITE "my.cdash.org")
  set(CTEST_DROP_LOCATION "/submit.php?project=Pele")
endif()

set(CTEST_DROP_SITE_CDASH TRUE)
