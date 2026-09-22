// A no-op falcON manipulator: request particle keys when reading the input.
// gyrfalcON's output "give=k" alone does not request keys during input.
#include <public/defman.h>

namespace falcON { namespace Manipulate {
class ocen_keep_keys : public manipulator {
public:
    ocen_keep_keys(const char*, const char*) {}
    const char* name() const { return "ocen_keep_keys"; }
    const char* describe() const { return "Preserve input particle keys; no dynamical action"; }
    fieldset need() const { return fieldset::k; }
    fieldset provide() const { return fieldset::empty; }
    fieldset change() const { return fieldset::empty; }
    bool manipulate(const snapshot*) const { return false; }
};
} }
__DEF__MAN__ALT(falcON::Manipulate::ocen_keep_keys);
