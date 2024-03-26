#include "utils.h"

void setTPSBlockList(char* p, int* countOut, size_t* TPSBlock) {
  int count = 0;
  while (*p != '\0') {
      // Skip any leading whitespace
      p += strspn(p, " \t");

      // Find the end of the number
      int len = strcspn(p, ", \t");

      // Temporarily replace the delimiter with '\0' to use atoi
      char temp = p[len];
      p[len] = '\0';

      // Convert the number and add it to the array
      TPSBlock[count++] = atoi(p);

      // Restore the original character
      p[len] = temp;

      // Skip to the next number
      p = strchr(p, ',');
      if (p == NULL) break; // No more numbers
      p++; // Skip the comma
  }
  *countOut = count;
}
