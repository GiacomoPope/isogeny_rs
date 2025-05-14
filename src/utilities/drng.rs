use rand_core::{CryptoRng, RngCore};
use sha2::{Digest, Sha512};

// Custom PRNG for tests.
// FOR TESTS ONLY. IT IS NOT NECESSARILY SECURE.
pub struct DRNG {
    buf: [u8; 64],
    ptr: usize,
}

impl DRNG {
    pub fn from_seed(seed: &[u8]) -> Self {
        let mut d = Self {
            buf: [0u8; 64],
            ptr: 0,
        };
        let mut sh = Sha512::new();
        sh.update(seed);
        d.buf[..].copy_from_slice(&sh.finalize());
        d
    }
}

impl RngCore for DRNG {
    fn next_u32(&mut self) -> u32 {
        let mut buf = [0u8; 4];
        self.fill_bytes(&mut buf);
        u32::from_le_bytes(buf)
    }

    fn next_u64(&mut self) -> u64 {
        let mut buf = [0u8; 8];
        self.fill_bytes(&mut buf);
        u64::from_le_bytes(buf)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        let len = dest.len();
        let mut off = 0;
        while off < len {
            let mut clen = 32 - self.ptr;
            if clen > (len - off) {
                clen = len - off;
            }
            dest[off..off + clen].copy_from_slice(&self.buf[self.ptr..self.ptr + clen]);
            self.ptr += clen;
            off += clen;
            if self.ptr == 32 {
                let mut sh = Sha512::new();
                sh.update(&self.buf);
                self.buf[..].copy_from_slice(&sh.finalize());
                self.ptr = 0;
            }
        }
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

impl CryptoRng for DRNG {}
