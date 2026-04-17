
#ifndef sph_test_h_included
#define sph_test_h_included

#include <stdio.h>
#include <sph-ga/sph/status.h>

#define test_helper_test(func) \
  printf("%s\n", #func); \
  status_require((func()))
#define test_helper_assert(description, expression) \
  if (!expression) { \
    printf("%s failed\n", description); \
    status_set_goto("", 1); \
  }
#define test_helper_assert_status(description, a) \
  status = a; \
  if (status_is_failure) { \
    printf("%s failed\n", description); \
    status_goto; \
  }
#define test_helper_display_summary_description(get_status_description) \
  if (status_is_success) { \
    printf(("--\ntests finished successfully.\n")); \
  } else { \
    printf(("\ntests failed. %d %s\n"), (status.id), (get_status_description(status))); \
  }
#define test_helper_display_summary \
  if (status_is_success) { \
    printf(("--\ntests finished successfully.\n")); \
  } else { \
    printf(("\ntests failed. %d\n"), (status.id)); \
  }

#define test_helper_assert1(description, expression, format_string, log_value) \
  if (!expression) { \
    printf(("%s failed. %s: " format_string "\n"), description, #log_value, log_value); \
    status_set_goto("", 1); \
  }
#define test_helper_assert2(description, expression, format_string1, log_value1, format_string2, log_value2) \
  if (!expression) { \
    printf(("%s failed. %s: " format_string1 ", %s: " format_string2 "\n"), description, #log_value1, log_value1, #log_value2, log_value2); \
    status_set_goto("", 1); \
  }
#define test_helper_assert3(description, expression, format_string1, log_value1, format_string2, log_value2, format_string3, log_value3) \
  if (!expression) { \
    printf(("%s failed. %s: "), description, #log_value1); \
    printf(format_string1, log_value1); \
    printf(", %s: ", #log_value2); \
    printf(format_string2, log_value2); \
    printf(", %s: ", #log_value3); \
    printf(format_string3, log_value3); \
    printf("\n"); \
    status_set_goto("", 1); \
  }
#define test_helper_assert4(description, expression, format_string1, log_value1, format_string2, log_value2, format_string3, log_value3, format_string4, log_value4) \
  if (!expression) { \
    printf(("%s failed. %s: "), description, #log_value1); \
    printf(format_string1, log_value1); \
    printf(", %s: ", #log_value2); \
    printf(format_string2, log_value2); \
    printf(", %s: ", #log_value3); \
    printf(format_string3, log_value3); \
    printf(", %s: ", #log_value4); \
    printf(format_string4, log_value4); \
    printf("\n"); \
    status_set_goto("", 1); \
  }
#endif
